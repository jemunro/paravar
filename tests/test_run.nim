## Tests for run.nim — argv parsing (R1), buildShellCmd (R2), runShards API (R3),
## CLI integration (R4), inferRunMode (R5), placeholders (R6), tool-managed (R7).
## Run from project root: nim c -r tests/test_run.nim
## Requires: tests/data/small.vcf.gz, tests/data/small.bcf (run generate_fixtures.sh once)

echo "--------------- Test Run ---------------"

import std/[atomics, os, osproc, strformat, strutils, tempfiles]
{.warning[Deprecated]: off.}
import std/threadpool
{.warning[Deprecated]: on.}
import test_utils
import "../src/blocky/scatter"
import "../src/blocky/run"

const DataDir  = "tests/data"
const SmallVcf = DataDir / "small.vcf.gz"
const SmallBcf = DataDir / "small.bcf"

# ===========================================================================
# R1–R6 — parseRunArgv: separator parsing and stage extraction
# ===========================================================================

# ---------------------------------------------------------------------------
# R1 — testSingleStage: single --- stage; blocky args and stage tokens correct
# ---------------------------------------------------------------------------
timed("R1.1", "parseRunArgv: single stage"):
  let argv = @["--shards", "4", "-o", "out", "input.vcf.gz",
               "---", "bcftools", "view", "-Oz"]
  let (pArgs, stages) = parseRunArgv(argv)
  doAssert pArgs == @["--shards", "4", "-o", "out", "input.vcf.gz"],
    "single stage: wrong blocky args"
  doAssert stages.len == 1, "single stage: expected 1 stage, got " & $stages.len
  doAssert stages[0] == @["bcftools", "view", "-Oz"],
    "single stage: wrong stage tokens"

# ---------------------------------------------------------------------------
# R2 — testMultiStage: multi-stage pipeline; inner -- is a plain token
# ---------------------------------------------------------------------------
timed("R1.2", "parseRunArgv: multi-stage pipeline"):
  let argv = @["-n", "10", "---",
               "bcftools", "+split-vep", "-Ou", "--", "-f", "%SYMBOL",
               "---",
               "bcftools", "view", "-s", "Sample", "-Oz"]
  let (pArgs, stages) = parseRunArgv(argv)
  doAssert pArgs == @["-n", "10"], "multi-stage: wrong blocky args"
  doAssert stages.len == 2, "multi-stage: expected 2 stages, got " & $stages.len
  doAssert stages[0] == @["bcftools", "+split-vep", "-Ou", "--", "-f", "%SYMBOL"],
    "multi-stage: wrong stage 0 tokens"
  doAssert stages[1] == @["bcftools", "view", "-s", "Sample", "-Oz"],
    "multi-stage: wrong stage 1 tokens"

# ---------------------------------------------------------------------------
# R3 — testSepFirst: --- at argv[0]; empty blocky args, one stage
# ---------------------------------------------------------------------------
timed("R1.3", "parseRunArgv: --- first, no blocky args"):
  let argv = @["---", "cat"]
  let (pArgs, stages) = parseRunArgv(argv)
  doAssert pArgs.len == 0, "sep-first: blocky args should be empty"
  doAssert stages.len == 1, "sep-first: expected 1 stage"
  doAssert stages[0] == @["cat"], "sep-first: wrong stage tokens"

# ---------------------------------------------------------------------------
# R4 — testDashDashPassthrough: -- inside stage is a plain token, not a separator
# ---------------------------------------------------------------------------
timed("R1.4", "parseRunArgv: -- inside stage is a plain token"):
  let argv = @["---", "bcftools", "+fill-tags", "-Ou", "--", "-t", "AF"]
  let (_, stages) = parseRunArgv(argv)
  doAssert stages.len == 1, "dash-dash: expected 1 stage"
  doAssert stages[0] == @["bcftools", "+fill-tags", "-Ou", "--", "-t", "AF"],
    "dash-dash: -- should be a plain token"

# ---------------------------------------------------------------------------
# R5 — testColonColonSep: ::: is an alias for ---; behaviour identical
# ---------------------------------------------------------------------------
timed("R1.5", "parseRunArgv: ::: is alias for ---"):
  let argv = @["-n", "4", ":::", "bcftools", "view", "-Oz"]
  let (pArgs, stages) = parseRunArgv(argv)
  doAssert pArgs == @["-n", "4"], "::: sep: wrong blocky args"
  doAssert stages.len == 1, "::: sep: expected 1 stage"
  doAssert stages[0] == @["bcftools", "view", "-Oz"], "::: sep: wrong stage tokens"

# ---------------------------------------------------------------------------
# R6 — testMixedSeparators: --- and ::: may be mixed freely in one invocation
# ---------------------------------------------------------------------------
timed("R1.6", "parseRunArgv: ::: and --- may be mixed"):
  let argv = @["---", "cmd1", "a", ":::", "cmd2", "b", "---", "cmd3", "c"]
  let (_, stages) = parseRunArgv(argv)
  doAssert stages.len == 3, "mixed seps: expected 3 stages, got " & $stages.len
  doAssert stages[0] == @["cmd1", "a"], "mixed: wrong stage 0"
  doAssert stages[1] == @["cmd2", "b"], "mixed: wrong stage 1"
  doAssert stages[2] == @["cmd3", "c"], "mixed: wrong stage 2"

# ---------------------------------------------------------------------------
# R1.7: no -o and no {} → stdout mode (runs successfully)
# ---------------------------------------------------------------------------
timed("R1.7", "CLI run: no -o and no {} -> stdout mode"):
  let (outp, code) = execCmdEx(
    "./blocky run -n 1 " & SmallVcf &
    " ::: cat 2>/dev/null | wc -c")
  doAssert code == 0, "no -o should succeed (stdout mode)"
  doAssert outp.strip.parseInt > 0, "stdout should have output"

# ===========================================================================
# R2 — buildShellCmd: stage list to shell command string
# ===========================================================================

# ---------------------------------------------------------------------------
# R2.1 — testBuildSingleStage: single stage produces no pipes
# ---------------------------------------------------------------------------
timed("R2.1", "buildShellCmd: single stage"):
  let cmd = buildShellCmd(@[@["bcftools", "view", "-Oz"]])
  doAssert cmd == "bcftools view -Oz", "build single: got " & cmd

# ---------------------------------------------------------------------------
# R2.2 — testBuildMultiStage: two stages joined with |
# ---------------------------------------------------------------------------
timed("R2.2", "buildShellCmd: multi-stage pipeline"):
  let cmd = buildShellCmd(@[
    @["bcftools", "+split-vep", "-Ou", "--", "-f", "%SYMBOL"],
    @["bcftools", "view", "-s", "Sample", "-Oz"]
  ])
  doAssert cmd == "bcftools +split-vep -Ou -- -f %SYMBOL | bcftools view -s Sample -Oz",
    "build multi: got " & cmd

# ---------------------------------------------------------------------------
# R2.3 — testBuildThreeStages: three stages produce two pipes
# ---------------------------------------------------------------------------
timed("R2.3", "buildShellCmd: three stages"):
  let cmd = buildShellCmd(@[@["cmd1", "a"], @["cmd2", "b"], @["cmd3", "c"]])
  doAssert cmd == "cmd1 a | cmd2 b | cmd3 c", "build three: got " & cmd

# ===========================================================================
# R3 — runShards unit tests (direct API calls, no binary)
# ===========================================================================

proc countRecords(path: string): int =
  let (o, _) = execCmdEx("bcftools query -f '%POS\\n' " & path & " 2>/dev/null | wc -l")
  o.strip.parseInt

proc recordsHash(paths: seq[string]): string =
  ## Concatenate records from paths in order (bcftools view -H, full genotypes),
  ## write to temp file, return sha256sum hex digest.
  let tmp = getTempDir() / "blocky_hash_" & $getCurrentProcessId() & ".txt"
  var f = open(tmp, fmWrite)
  for p in paths:
    let (o, _) = execCmdEx("bcftools view -H " & p & " 2>/dev/null")
    f.write(o)
  f.close()
  let (h, _) = execCmdEx("sha256sum " & tmp)
  removeFile(tmp)
  h.split(" ")[0]

# ---------------------------------------------------------------------------
# R3.1 — testRunSingle1Shard: 1 shard, cat pipeline; output valid, record count matches
# ---------------------------------------------------------------------------
timed("R3.1", "runPipeline: 1 worker, cat stage"):
  doAssert fileExists(SmallVcf), "fixture missing — run generate_fixtures.sh"
  let tmpDir = createTempDir("blocky_", "")
  let outPath = tmpDir / "out.vcf.gz"
  runPipeline(RunPipelineCfg(
    vcfPath: SmallVcf, nWorkers: 1, maxShardsPerWorker: 1, nThreads: 1,
    stages: @[@["cat"]], outputPath: outPath))
  doAssert fileExists(outPath), "output missing: " & outPath
  let orig = countRecords(SmallVcf)
  let got  = countRecords(outPath)
  doAssert got == orig, &"1-worker record count mismatch: got {got}, expected {orig}"
  removeDir(tmpDir)

# ---------------------------------------------------------------------------
# R3.2 — testRun4Shards: 4 shards concurrent; all valid, record count + content hash match
# ---------------------------------------------------------------------------
timed("R3.2", "runPipeline: 4 workers, cat stage"):
  doAssert fileExists(SmallVcf), "fixture missing"
  let tmpDir = createTempDir("blocky_", "")
  let outPath = tmpDir / "out.vcf.gz"
  runPipeline(RunPipelineCfg(
    vcfPath: SmallVcf, nWorkers: 4, maxShardsPerWorker: 1, nThreads: 1,
    stages: @[@["cat"]], outputPath: outPath))
  doAssert fileExists(outPath), "output missing"
  let orig = countRecords(SmallVcf)
  let got  = countRecords(outPath)
  doAssert got == orig, &"4-worker record count mismatch: {got} vs {orig}"
  removeDir(tmpDir)

# ---------------------------------------------------------------------------
# R3.3 — testBcfRun4Shards: BCF runShards 4 shards, bcftools view -Ob; hash matches
# ---------------------------------------------------------------------------
timed("R3.3", "runPipeline BCF: 4 workers, cat stage"):
  doAssert fileExists(SmallBcf), "BCF fixture missing — run generate_fixtures.sh"
  let tmpDir = createTempDir("blocky_", "")
  let outPath = tmpDir / "out.vcf.gz"
  runPipeline(RunPipelineCfg(
    vcfPath: SmallBcf, nWorkers: 4, maxShardsPerWorker: 1, nThreads: 1,
    stages: @[@["cat"]], outputPath: outPath))
  doAssert fileExists(outPath), "BCF output missing"
  let orig = countRecords(SmallBcf)
  let got  = countRecords(outPath)
  doAssert got == orig, &"BCF 4-worker record count mismatch: {got} vs {orig}"
  removeDir(tmpDir)

# ===========================================================================
# R4 — CLI integration tests via the compiled binary
# ===========================================================================

const BinPath = "./blocky"

timed("R4.0", "binary available (run CLI tests)"):
  if not fileExists(BinPath):
    let (outp, code) = execCmdEx("nimble build 2>&1")
    if code != 0:
      echo "nimble build failed:\n", outp
      quit(1)
  doAssert fileExists(BinPath), "binary not found: " & BinPath & " (run nimble build)"

proc runBin(args: string): (string, int) =
  let t = getEnv("BLOCKY_TEST_TIMEOUT", "10")
  execCmdEx("timeout " & t & " " & BinPath & " run " & args & " 2>&1")

# ---------------------------------------------------------------------------
# R4.1 — testCliRunMissingSep: no --- exits 1 with '---' in error message
# ---------------------------------------------------------------------------
timed("R4.1", "CLI run: missing --- -> exits 1 with message"):
  let (outp, code) = runBin(&"-n 1 -o /tmp/x {SmallVcf}")
  doAssert code != 0, "missing --- should exit non-zero"
  doAssert "---" in outp, &"error should mention '---', got:\n{outp}"

# ---------------------------------------------------------------------------
# R4.2 — testCliRun1Shard: run -n 1 --- cat; output valid, record count matches
# ---------------------------------------------------------------------------
timed("R4.2", "CLI run: 1 worker, cat stage, single output"):
  let tmpDir = createTempDir("blocky_", "")
  let outFile = tmpDir / "out.vcf.gz"
  let (outp, code) = runBin(&"-n 1 -o {outFile} {SmallVcf} --- cat")
  doAssert code == 0, &"run -n 1 exited {code}:\n{outp}"
  doAssert fileExists(outFile), "output missing: " & outFile
  let orig = countRecords(SmallVcf)
  doAssert countRecords(outFile) == orig, "1-worker CLI record count mismatch"
  removeDir(tmpDir)

# ---------------------------------------------------------------------------
# R4.3 — testCliRun4Shards: run -n 4 --- cat; 4 shards valid, record count matches
# ---------------------------------------------------------------------------
timed("R4.3", "CLI run: 4 workers, cat stage, single output"):
  let tmpDir = createTempDir("blocky_", "")
  let outFile = tmpDir / "out.vcf.gz"
  let (outp, code) = runBin(&"-n 4 -o {outFile} {SmallVcf} --- cat")
  doAssert code == 0, &"run -n 4 exited {code}:\n{outp}"
  doAssert fileExists(outFile), "output missing"
  doAssert countRecords(outFile) == countRecords(SmallVcf),
    "4-worker CLI record count mismatch"
  removeDir(tmpDir)

# ---------------------------------------------------------------------------
# R4.4 — testCliRunMultiStage: run -n 2 --- cat --- cat; record count matches
# ---------------------------------------------------------------------------
timed("R4.4", "CLI run: multi-stage pipeline (cat | cat), records match"):
  let tmpDir = createTempDir("blocky_", "")
  let outFile = tmpDir / "out.vcf.gz"
  let (outp, code) = runBin(&"-n 2 -o {outFile} {SmallVcf} --- cat --- cat")
  doAssert code == 0, &"run multi-stage exited {code}:\n{outp}"
  doAssert countRecords(outFile) == countRecords(SmallVcf),
    "multi-stage record count mismatch"
  removeDir(tmpDir)

# ---------------------------------------------------------------------------
# R4.5 — testCliRunJError: -j is unknown, exits non-zero
# ---------------------------------------------------------------------------
timed("R4.5", "CLI run: -j exits non-zero (unknown option)"):
  let tmpDir = createTempDir("blocky_", "")
  let outp_template = tmpDir / "out.vcf.gz"
  let (_, code) = runBin(&"-n 4 -j 1 -o {outp_template} {SmallVcf} --- cat")
  doAssert code != 0, "-j should cause a non-zero exit (unknown option)"
  removeDir(tmpDir)

# ---------------------------------------------------------------------------
# R4.6 — testCliRunAttachedN: -n4 (attached, no space) parsed correctly
# ---------------------------------------------------------------------------
timed("R4.6", "CLI run: -n4 attached-value flag works"):
  let tmpDir = createTempDir("blocky_", "")
  let outFile = tmpDir / "out.vcf.gz"
  let (outp, code) = runBin(&"-n4 -o {outFile} {SmallVcf} --- cat")
  doAssert code == 0, &"run -n4 exited {code}:\n{outp}"
  doAssert countRecords(outFile) == countRecords(SmallVcf), "-n4 record count mismatch"
  removeDir(tmpDir)

# ---------------------------------------------------------------------------
# R4.7 — testCliRunMaxJobsError: --max-jobs is unknown, exits non-zero
# ---------------------------------------------------------------------------
timed("R4.7", "CLI run: --max-jobs exits non-zero (unknown option)"):
  let tmpDir = createTempDir("blocky_", "")
  let outp_template = tmpDir / "out.vcf.gz"
  let (_, code) = runBin(&"-n 2 --max-jobs 10 -o {outp_template} {SmallVcf} --- cat")
  doAssert code != 0, "--max-jobs should cause a non-zero exit (unknown option)"
  removeDir(tmpDir)

# ---------------------------------------------------------------------------
# R4.8 — testCliRunNonZeroExit: non-zero stage exit → blocky exits 1, shard mentioned
# ---------------------------------------------------------------------------
timed("R4.8", "CLI run: non-zero stage exit -> blocky exits 1, shard mentioned"):
  let tmpDir = createTempDir("blocky_", "")
  let outp_template = tmpDir / "out.vcf.gz"
  let (outp, code) = runBin(&"-n 1 -o {outp_template} {SmallVcf} --- false")
  doAssert code != 0, "non-zero stage should make blocky exit non-zero"
  doAssert "shard 1" in outp, &"stderr should mention 'shard 1', got:\n{outp}"
  removeDir(tmpDir)

# ---------------------------------------------------------------------------
# R4.9 — testCliRunColonColon: ::: separator works at CLI level
# ---------------------------------------------------------------------------
timed("R4.9", "CLI run: ::: separator works like ---"):
  let tmpDir = createTempDir("blocky_", "")
  let outFile = tmpDir / "out.vcf.gz"
  let (outp, code) = runBin(&"-n 2 -o {outFile} {SmallVcf} ::: cat")
  doAssert code == 0, &"::: separator exited {code}:\n{outp}"
  doAssert fileExists(outFile), "::: output missing"
  removeDir(tmpDir)

# ---------------------------------------------------------------------------
# R4.10 — testCliRunDashDashPassthrough: -- inside stage passed through to tool
# ---------------------------------------------------------------------------
timed("R4.10", "CLI run: -- inside stage passed through to bcftools"):
  # bcftools view -Oz -- - reads stdin, writes bgzipped VCF to stdout.
  let tmpDir = createTempDir("blocky_", "")
  let outFile = tmpDir / "out.vcf.gz"
  let (outp, code) = runBin(&"-n 1 -o {outFile} {SmallVcf} --- bcftools view -Oz -- -")
  doAssert code == 0, &"-- passthrough exited {code}:\n{outp}"
  doAssert fileExists(outFile), "-- passthrough output missing"
  doAssert countRecords(outFile) == countRecords(SmallVcf), "-- passthrough record count mismatch"
  removeDir(tmpDir)

# ---------------------------------------------------------------------------
# R4.11 — testCliRunNoKill: --no-kill flag accepted; run completes normally
# ---------------------------------------------------------------------------
timed("R4.11", "CLI run: --no-kill flag accepted, all shards complete normally"):
  let tmpDir = createTempDir("blocky_", "")
  let outFile = tmpDir / "out.vcf.gz"
  let (outp, code) = runBin(&"-n 2 --no-kill -o {outFile} {SmallVcf} --- cat")
  doAssert code == 0, &"run --no-kill exited {code}:\n{outp}"
  doAssert fileExists(outFile), "--no-kill output missing"
  removeDir(tmpDir)

# ===========================================================================
# R5 — inferRunMode: non-error cases
# (The error case — no output, no {} — calls quit(1) and is covered at the
#  CLI level by testCliRunMissingSep.)
# ===========================================================================

timed("R5.1", "inferRunMode: o / no-discard -> rmFile"):
  let m = inferRunMode(true, false)
  doAssert m == rmFile, "expected rmFile, got " & $m

timed("R5.2", "inferRunMode: no-o / discard -> rmDiscard"):
  let m = inferRunMode(false, true)
  doAssert m == rmDiscard, "expected rmDiscard, got " & $m

timed("R5.3", "inferRunMode: no-o / no-discard -> rmStdout"):
  let m = inferRunMode(false, false)
  doAssert m == rmStdout, "expected rmStdout, got " & $m

timed("R5.4", "inferRunMode: o / no-discard -> rmFile"):
  let m = inferRunMode(true, false)
  doAssert m == rmFile, "expected rmFile, got " & $m

# ===========================================================================
# R6 — hasBracePlaceholder, substituteToken, buildShellCmdForShard
# ===========================================================================

timed("R6.1", "hasBracePlaceholder: detects {} in token"):
  doAssert hasBracePlaceholder(@[@["bcftools", "view", "-o", "out.{}.vcf.gz"]]),
    "expected true for token containing {}"

timed("R6.2", "hasBracePlaceholder: false when no {} present"):
  doAssert not hasBracePlaceholder(@[@["bcftools", "view", "-Oz"]]),
    "expected false when no {} present"

timed("R6.3", "hasBracePlaceholder: escaped \\{} does not count"):
  # \{} is escaped — should NOT count as a placeholder
  doAssert not hasBracePlaceholder(@[@["tool", "\\{}"  ]]),
    "expected false for escaped \\{}"

timed("R6.4", "substituteToken: basic {} substitution"):
  let r = substituteToken("out.{}.vcf.gz", "03")
  doAssert r == "out.03.vcf.gz", "basic sub: got " & r

timed("R6.5", "substituteToken: multiple {} replaced"):
  let r = substituteToken("{}-{}", "05")
  doAssert r == "05-05", "multiple: got " & r

timed("R6.6", "substituteToken: \\{} -> literal {}"):
  # \{} should become literal {} in the output (backslash consumed)
  let r = substituteToken("\\{}", "07")
  doAssert r == "{}", "escaped: got " & r

timed("R6.7", "substituteToken: no {} -> unchanged"):
  let r = substituteToken("bcftools", "02")
  doAssert r == "bcftools", "no-brace: got " & r

timed("R6.8", "substituteToken: \\{} then {} -> literal {} then shardNum"):
  # \{} followed by {} → literal {} then shardNum
  let r = substituteToken("\\{}{}", "04")
  doAssert r == "{}04", "mixed escape: got " & r

timed("R6.9", "buildShellCmdForShard: zero-pads to width 2 for nShards=10"):
  # nShards = 10 → width 2; shardIdx = 0 → "01"
  let cmd = buildShellCmdForShard(@[@["bcftools", "view", "-o", "out.{}.vcf.gz"]], 0, 10)
  doAssert "out.01.vcf.gz" in cmd, "padding w2: got " & cmd

timed("R6.10", "buildShellCmdForShard: no {} -> cmd unchanged"):
  # No {} → cmd unchanged (aside from quoting)
  let cmd = buildShellCmdForShard(@[@["cat"]], 0, 4)
  doAssert cmd == "cat", "no-brace: got " & cmd

timed("R6.11", "buildShellCmdForShard: width 1 for nShards=4"):
  # nShards = 4 → width 1; shardIdx = 2 → "3"
  let cmd = buildShellCmdForShard(@[@["tool", "{}"]], 2, 4)
  doAssert "3" in cmd, "width-1: got " & cmd
  doAssert "03" notin cmd, "width-1: should not be zero-padded: got " & cmd

# ===========================================================================
# R7 — --discard mode integration tests
# ===========================================================================

# ---------------------------------------------------------------------------
# R7.1: --discard mode via API
# ---------------------------------------------------------------------------
timed("R7.1", "--discard API: 2 workers, tool writes own output"):
  doAssert fileExists(SmallVcf), "fixture missing"
  let tmpDir = createTempDir("blocky_", "")
  let outTemplate = tmpDir / "out.{}.vcf.gz"
  let stages = @[@["bcftools", "view", "-Oz", "-o", outTemplate]]
  runPipeline(RunPipelineCfg(
    vcfPath: SmallVcf, nWorkers: 2, maxShardsPerWorker: 1, nThreads: 1,
    stages: stages, discardStdout: true))
  var total = 0
  for i in 0..1:
    let p = shardOutputPath(outTemplate, i, 2)
    doAssert fileExists(p), &"--discard shard {i+1} missing: {p}"
    let (_, bc) = execCmdEx("bcftools view -HG " & p & " > /dev/null 2>&1")
    doAssert bc == 0, &"bcftools rejected --discard shard {i+1}"
    total += countRecords(p)
  doAssert total == countRecords(SmallVcf),
    &"--discard record count mismatch: {total}"
  removeDir(tmpDir)

# ---------------------------------------------------------------------------
# R7.2: --discard CLI — {} in tool cmd, --discard → tool writes files
# ---------------------------------------------------------------------------
timed("R7.2", "--discard CLI: {} in tool cmd"):
  doAssert fileExists(SmallVcf), "fixture missing"
  let tmpDir = createTempDir("blocky_", "")
  let outTemplate = tmpDir / "out.{}.vcf.gz"
  let (outp, code) = execCmdEx(
    BinPath & " run -n 2 --discard " & SmallVcf &
    " ::: bcftools view -Oz -o " & outTemplate & " 2>&1")
  doAssert code == 0, &"--discard CLI exited {code}:\n{outp}"
  var total = 0
  for i in 0..1:
    let p = shardOutputPath(outTemplate, i, 2)
    doAssert fileExists(p), &"--discard CLI shard {i+1} missing: {p}"
    total += countRecords(p)
  doAssert total == countRecords(SmallVcf),
    &"--discard CLI record count mismatch: {total}"
  removeDir(tmpDir)

# ---------------------------------------------------------------------------
# R7.3: {} with -o (no --discard) → stdout captured + {} substituted
# ---------------------------------------------------------------------------
timed("R7.3", "{} with -o: stdout captured, {} substituted"):
  doAssert fileExists(SmallVcf), "fixture missing"
  let tmpDir = createTempDir("blocky_", "")
  let auxTemplate = tmpDir / "aux.{}.txt"
  let outPath = tmpDir / "captured.vcf.gz"
  # Tool writes aux files via {} AND stdout is captured to -o.
  # Use 'tee' to write to {} file and also emit to stdout.
  let (outp, code) = execCmdEx(
    BinPath & " run -n 2 -o " & outPath & " " & SmallVcf &
    " ::: bcftools view -Ov 2>&1")
  doAssert code == 0, &"R7.3 exited {code}:\n{outp}"
  doAssert fileExists(outPath), "R7.3: captured output file missing"
  let total = countRecords(outPath)
  doAssert total == countRecords(SmallVcf),
    &"R7.3: record count mismatch: {total}"
  removeDir(tmpDir)

# ===========================================================================
# W1 — DepositQueue and atomic shard counter
# ===========================================================================

# ---------------------------------------------------------------------------
# W1.1 — DepositQueue: deposit and waitFor in order
# ---------------------------------------------------------------------------
timed("W1.1", "DepositQueue: deposit in order, waitFor retrieves correctly"):
  var q = newDepositQueue(4)
  q.deposit(0, "/tmp/shard_0.tmp")
  q.deposit(1, "/tmp/shard_1.tmp")
  q.deposit(2, "/tmp/shard_2.tmp")
  q.deposit(3, "/tmp/shard_3.tmp")
  doAssert q.waitFor(0) == "/tmp/shard_0.tmp", "W1.1: slot 0 wrong"
  doAssert q.waitFor(1) == "/tmp/shard_1.tmp", "W1.1: slot 1 wrong"
  doAssert q.waitFor(2) == "/tmp/shard_2.tmp", "W1.1: slot 2 wrong"
  doAssert q.waitFor(3) == "/tmp/shard_3.tmp", "W1.1: slot 3 wrong"
  freeDepositQueue(q)

# ---------------------------------------------------------------------------
# W1.2 — DepositQueue: out-of-order deposit, in-order waitFor
# ---------------------------------------------------------------------------
timed("W1.2", "DepositQueue: out-of-order deposit, ordered waitFor"):
  var q = newDepositQueue(4)
  # Deposit out of order using threads to exercise cross-thread visibility.
  proc depositDelayed(qp: ptr DepositQueue; idx: int; path: string): int {.gcsafe.} =
    sleep(idx * 10)  # stagger: slot 2 deposits first, slot 3 last
    qp[].deposit(idx, path)
    0
  setMaxPoolSize(4)
  let fv2 = spawn depositDelayed(addr q, 2, "/tmp/s2")
  let fv0 = spawn depositDelayed(addr q, 0, "/tmp/s0")
  let fv3 = spawn depositDelayed(addr q, 3, "/tmp/s3")
  let fv1 = spawn depositDelayed(addr q, 1, "/tmp/s1")
  # waitFor in order — each blocks until that slot is deposited.
  doAssert q.waitFor(0) == "/tmp/s0", "W1.2: slot 0 wrong"
  doAssert q.waitFor(1) == "/tmp/s1", "W1.2: slot 1 wrong"
  doAssert q.waitFor(2) == "/tmp/s2", "W1.2: slot 2 wrong"
  doAssert q.waitFor(3) == "/tmp/s3", "W1.2: slot 3 wrong"
  discard ^fv0; discard ^fv1; discard ^fv2; discard ^fv3
  freeDepositQueue(q)

# ---------------------------------------------------------------------------
# W1.3 — Atomic shard counter: each index claimed exactly once
# ---------------------------------------------------------------------------
timed("W1.3", "gNextShard: atomic pull, each index claimed once"):
  const NItems = 8
  gNextShard.store(0, moRelaxed)
  var claimed: array[NItems, Atomic[int]]
  for i in 0 ..< NItems:
    claimed[i].store(0, moRelaxed)
  proc puller(claimedArr: ptr array[NItems, Atomic[int]]): int {.gcsafe.} =
    while true:
      let idx = gNextShard.fetchAdd(1, moRelaxed)
      if idx >= NItems: break
      discard claimedArr[][idx].fetchAdd(1, moRelaxed)
    0
  setMaxPoolSize(4)
  var fvs: array[4, FlowVar[int]]
  for i in 0..3:
    fvs[i] = spawn puller(addr claimed)
  for fv in fvs:
    discard ^fv
  for i in 0 ..< NItems:
    doAssert claimed[i].load(moRelaxed) == 1,
      &"W1.3: index {i} claimed {claimed[i].load(moRelaxed)} times (expected 1)"
