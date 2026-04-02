## run — scatter a VCF into N shards and pipe each through a tool pipeline.
##
## This module is responsible for:
##   1. Parsing the "---"-separated argv into paravar args + pipeline stages.
##   2. Building the sh -c command string for each shard.
##   3. Executing per-shard pipelines concurrently up to --jobs.

import std/[cpuinfo, os, posix, sequtils, strformat, strutils]
{.warning[Deprecated]: off.}
import std/threadpool
{.warning[Deprecated]: on.}
import scatter
import gather

# ---------------------------------------------------------------------------
# Argv parsing
# ---------------------------------------------------------------------------

proc isSep(tok: string): bool {.inline.} =
  ## Return true for "---" or ":::" — both are valid pipeline stage separators.
  tok == "---" or tok == ":::"

proc parseRunArgv*(argv: seq[string]): (seq[string], seq[seq[string]]) =
  ## Split argv at "---" or ":::" separators (both are equivalent).
  ## Returns (paravarArgs, stages) where stages is a seq of per-stage token seqs.
  ## Exits 1 with a message if no separator is present or if any stage is empty.
  var firstSep = -1
  for i, tok in argv:
    if isSep(tok):
      firstSep = i
      break
  if firstSep < 0:
    stderr.writeLine "paravar run: at least one --- stage is required"
    quit(1)
  let paravarArgs = argv[0 ..< firstSep]
  var stages: seq[seq[string]]
  var cur: seq[string]
  for i in firstSep + 1 ..< argv.len:
    if isSep(argv[i]):
      if cur.len == 0:
        stderr.writeLine "paravar run: empty pipeline stage"
        quit(1)
      stages.add(cur)
      cur = @[]
    else:
      cur.add(argv[i])
  if cur.len == 0:
    stderr.writeLine "paravar run: empty pipeline stage"
    quit(1)
  stages.add(cur)
  result = (paravarArgs, stages)

# ---------------------------------------------------------------------------
# Shell command construction
# ---------------------------------------------------------------------------

proc buildShellCmd*(stages: seq[seq[string]]): string =
  ## Build a sh -c command string from pipeline stages.
  ## Each token is shell-quoted so special characters (< > | & etc.) in
  ## filter expressions are passed through safely.  Stages are joined with " | ".
  var parts: seq[string]
  for stage in stages:
    parts.add(stage.mapIt(quoteShell(it)).join(" "))
  result = parts.join(" | ")

# ---------------------------------------------------------------------------
# Per-shard pipe execution (R3/R4)
# ---------------------------------------------------------------------------

type InFlight = object
  ## Tracks one active shard: child process + writer thread.
  pid:      Pid
  writeFv:  FlowVar[int]
  shardIdx: int

proc forkExecSh(pipeReadFd: cint; pipeWriteFd: cint; stdoutFd: cint;
                shellCmd: string; shardIdx: int): Pid =
  ## Fork a child that runs sh -c shellCmd with stdin = pipeReadFd and
  ## stdout = stdoutFd.  stderr is inherited.  Returns child PID.
  ## pipeWriteFd is closed in the child so the child does not hold the
  ## write-end of its own stdin pipe (which would prevent EOF).
  let pid = posix.fork()
  if pid < 0:
    stderr.writeLine &"error: fork() failed for shard {shardIdx + 1}"
    quit(1)
  if pid == 0:
    # Child: rewire stdin and stdout, close the pipe write-end, exec shell.
    if posix.dup2(pipeReadFd, STDIN_FILENO) < 0 or
       posix.dup2(stdoutFd,   STDOUT_FILENO) < 0:
      exitnow(1)
    discard posix.close(pipeReadFd)
    discard posix.close(pipeWriteFd)
    discard posix.close(stdoutFd)
    let args = allocCStringArray(["sh", "-c", shellCmd])
    discard posix.execvp("sh", args)
    deallocCStringArray(args)
    exitnow(127)
  result = pid

proc killAll(running: seq[InFlight]) =
  ## Send SIGTERM to every in-flight child process.
  for s in running:
    discard posix.kill(s.pid, SIGTERM)

proc waitOne(running: var seq[InFlight]; failed: var bool) =
  ## Wait for any one child to finish; record failure if exit code != 0.
  var status: cint
  let donePid = posix.waitpid(-1, status, 0)
  let code    = int((status shr 8) and 0xff)
  var j = 0
  while j < running.len:
    if running[j].pid == donePid:
      discard ^running[j].writeFv
      if code != 0:
        stderr.writeLine &"shard {running[j].shardIdx + 1}: pipeline exited with code {code}"
        failed = true
      running.del(j)
      return
    j += 1

proc runShards*(vcfPath: string; nShards: int; outputTemplate: string;
                nThreads: int; forceScan: bool; maxJobs: int;
                shellCmd: string; noKill: bool = false) =
  ## Scatter vcfPath into nShards shards and pipe each through shellCmd.
  ## outputTemplate may contain {} or use shard_NN. prefix naming.
  ## Up to maxJobs shards run concurrently; pass 0 to use all CPUs.
  ## On failure, siblings are killed (SIGTERM) unless noKill is true.
  ## If any shard pipeline exits non-zero, prints to stderr and exits 1 at end.
  let actualThreads = if nThreads == 0: countProcessors() else: nThreads
  let actualMaxJobs = if maxJobs == 0: countProcessors() else: maxJobs
  setMaxPoolSize(max(actualThreads, actualMaxJobs))
  let fmt      = if vcfPath.endsWith(".bcf"): FileFormat.Bcf else: FileFormat.Vcf
  let tasks    = computeShards(vcfPath, nShards, actualThreads, forceScan, fmt)
  var anyFailed = false
  var inFlight: seq[InFlight]
  for i in 0 ..< nShards:
    # Stop launching new shards if a failure has been detected and !noKill.
    if anyFailed and not noKill: break
    # Drain one finished shard if at job capacity.
    while inFlight.len >= actualMaxJobs:
      waitOne(inFlight, anyFailed)
      if anyFailed and not noKill: break
    if anyFailed and not noKill: break
    let outPath = shardOutputPath(outputTemplate, i, nShards)
    createDir(outPath.parentDir)
    # Create pipe: [0] = read-end (child stdin), [1] = write-end (shard writer).
    var pipeFds: array[2, cint]
    if posix.pipe(pipeFds) != 0:
      stderr.writeLine &"error: pipe() failed for shard {i + 1}"
      quit(1)
    let outFileFd = posix.open(outPath.cstring,
                               O_WRONLY or O_CREAT or O_TRUNC,
                               0o666.Mode)
    if outFileFd < 0:
      stderr.writeLine &"error: could not create output file: {outPath}"
      quit(1)
    let pid = forkExecSh(pipeFds[0], pipeFds[1], outFileFd, shellCmd, i)
    # Parent closes the fds the child owns.
    discard posix.close(pipeFds[0])
    discard posix.close(outFileFd)
    # Assign pipe write-end to the shard task and spawn the writer thread.
    var task = tasks[i]
    task.outFd = pipeFds[1]
    let writeFv = spawn doWriteShard(task)
    inFlight.add(InFlight(pid: pid, writeFv: writeFv, shardIdx: i))
  # On failure without --no-kill: terminate all remaining children, then reap.
  if anyFailed and not noKill:
    killAll(inFlight)
  while inFlight.len > 0:
    waitOne(inFlight, anyFailed)
  if anyFailed:
    quit(1)

# ---------------------------------------------------------------------------
# G6 — gather path: interceptor threads + concatenation
# ---------------------------------------------------------------------------

type InFlightGather = object
  ## Tracks one gather shard: shell process + writer thread + interceptor thread.
  pid:         Pid
  writeFv:     FlowVar[int]
  interceptFv: FlowVar[int]
  shardIdx:    int
  tmpPath:     string

proc killAllGather(running: seq[InFlightGather]) =
  for s in running:
    discard posix.kill(s.pid, SIGTERM)

proc waitOneGather(running: var seq[InFlightGather]; failed: var bool) =
  ## Wait for any one shell child to exit; sync its writer and interceptor; record failure.
  var status: cint
  let donePid = posix.waitpid(-1, status, 0)
  let code    = int((status shr 8) and 0xff)
  var j = 0
  while j < running.len:
    if running[j].pid == donePid:
      discard ^running[j].writeFv
      let interceptCode = ^running[j].interceptFv
      if code != 0 or interceptCode != 0:
        stderr.writeLine &"shard {running[j].shardIdx + 1}: pipeline exited with code {code}"
        failed = true
      running.del(j)
      return
    j += 1

proc runShardsGather*(vcfPath: string; nShards: int; outputTemplate: string;
                      nThreads: int; forceScan: bool; maxJobs: int;
                      shellCmd: string; noKill: bool; cfg: GatherConfig) =
  ## Scatter vcfPath into nShards shards, pipe each through shellCmd, capture stdout
  ## via interceptor threads, then concatenate temp files into cfg.outputPath.
  let actualThreads = if nThreads == 0: countProcessors() else: nThreads
  let actualMaxJobs = if maxJobs == 0: countProcessors() else: maxJobs
  # Pool must accommodate scatter writer threads and interceptor threads simultaneously.
  setMaxPoolSize(max(actualThreads, actualMaxJobs) + actualMaxJobs)
  let fmt   = if vcfPath.endsWith(".bcf"): FileFormat.Bcf else: FileFormat.Vcf
  let tasks = computeShards(vcfPath, nShards, actualThreads, forceScan, fmt)
  createDir(cfg.tmpDir)
  var anyFailed = false
  var inFlight:    seq[InFlightGather]
  var allTmpPaths: seq[string]
  for i in 0 ..< nShards:
    if anyFailed and not noKill: break
    while inFlight.len >= actualMaxJobs:
      waitOneGather(inFlight, anyFailed)
      if anyFailed and not noKill: break
    if anyFailed and not noKill: break
    let tmpPath =
      if i == 0: cfg.outputPath
      else:
        let shardBase = shardOutputPath(cfg.outputPath, i, nShards).lastPathPart
        cfg.tmpDir / "paravar_" & shardBase & ".tmp"
    if i > 0:
      allTmpPaths.add(tmpPath)
    # stdin pipe:  shard writer → shell stdin
    # stdout pipe: shell stdout → interceptor
    var stdinPipe:  array[2, cint]
    var stdoutPipe: array[2, cint]
    if posix.pipe(stdinPipe) != 0 or posix.pipe(stdoutPipe) != 0:
      stderr.writeLine &"error: pipe() failed for shard {i + 1}"
      quit(1)
    let pid = forkExecSh(stdinPipe[0], stdinPipe[1], stdoutPipe[1], shellCmd, i)
    # Parent closes fds that now belong to the child.
    discard posix.close(stdinPipe[0])
    discard posix.close(stdoutPipe[1])
    # Spawn shard writer (writes shard bytes to stdinPipe[1]).
    var task = tasks[i]
    task.outFd = stdinPipe[1]
    let writeFv = spawn doWriteShard(task)
    # Spawn interceptor (reads shell stdout from stdoutPipe[0], writes to tmpPath).
    let interceptFd = stdoutPipe[0]
    var cfgCopy = cfg
    let interceptFv = spawn runInterceptor(cfgCopy, i, interceptFd, tmpPath)
    inFlight.add(InFlightGather(pid: pid, writeFv: writeFv, interceptFv: interceptFv,
                                shardIdx: i, tmpPath: tmpPath))
  if anyFailed and not noKill:
    killAllGather(inFlight)
  while inFlight.len > 0:
    waitOneGather(inFlight, anyFailed)
  if anyFailed:
    cleanupTempDir(cfg.tmpDir, allTmpPaths, false)
    quit(1)
  concatenateShards(cfg, allTmpPaths)
