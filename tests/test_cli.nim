## Tests for CLI argument handling.
## Run from project root: nim c -r tests/test_cli.nim
## Requires the vcfparty binary to be built first (nimble build).

echo "--------------- Test CLI ---------------"

import std/[os, osproc, sequtils, strformat, strutils, tempfiles]
import test_utils

const BinPath  = "./vcfparty"
const DataDir  = "tests/data"
const SmallVcf = DataDir / "small.vcf.gz"       # TBI indexed
const CsiVcf   = DataDir / "small_csi.vcf.gz"   # CSI indexed only
const KgVcf    = DataDir / "chr22_1kg.vcf.gz"
const SmallBcf = DataDir / "small.bcf"          # CSI indexed BCF

proc run(args: string): (string, int) =
  ## Run partyvcf with shell args; combine stdout+stderr; return (outp, code).
  execCmdEx(BinPath & " " & args & " 2>&1")

proc recordsHash(paths: seq[string]): string =
  ## Concatenate records from paths in order (bcftools view -H, full genotypes),
  ## write to temp file, return sha256sum hex digest.
  let tmp = getTempDir() / "vcfparty_hash_" & $getCurrentProcessId() & ".txt"
  var f = open(tmp, fmWrite)
  for p in paths:
    let (o, _) = execCmdEx("bcftools view -H " & p & " 2>/dev/null")
    f.write(o)
  f.close()
  let (h, _) = execCmdEx("sha256sum " & tmp)
  removeFile(tmp)
  h.split(" ")[0]

# ---------------------------------------------------------------------------
# Build binary (setup)
# ---------------------------------------------------------------------------

timed("CL0", "binary available"):
  if not fileExists(BinPath):
    let (outp, code) = execCmdEx("nimble build 2>&1")
    if code != 0:
      echo "nimble build failed:\n", outp
      quit(1)
  doAssert fileExists(BinPath), "binary not found: " & BinPath & " (run nimble build)"

# ===========================================================================
# CL1–CL7 — Error cases (missing flags, invalid args, unsupported modes)
# ===========================================================================

# ---------------------------------------------------------------------------
# CL1 — testMissingN: missing -n exits non-zero with 'n' in error message
# ---------------------------------------------------------------------------

timed("CL1", "missing -n exits non-zero"):
  let (outp, code) = run(&"scatter -o /tmp/vcfparty_cli_test {SmallVcf}")
  doAssert code != 0, "missing -n should exit non-zero"
  doAssert "n" in outp.toLowerAscii,
    &"missing -n error should mention 'n', got: {outp}"

# ---------------------------------------------------------------------------
# CL2 — testInvalidN0: -n 0 exits non-zero
# ---------------------------------------------------------------------------

timed("CL2", "-n 0 exits non-zero"):
  let (_, code) = run(&"scatter -n 0 -o /tmp/vcfparty_cli_test {SmallVcf}")
  doAssert code != 0, "-n 0 should exit non-zero"

# ---------------------------------------------------------------------------
# CL3 — testMissingO: missing -o exits non-zero with 'o' in error message
# ---------------------------------------------------------------------------

timed("CL3", "missing -o exits non-zero"):
  let (outp, code) = run(&"scatter -n 2 {SmallVcf}")
  doAssert code != 0, "missing -o should exit non-zero"
  doAssert "o" in outp.toLowerAscii,
    &"missing -o error should mention 'o', got: {outp}"

# ---------------------------------------------------------------------------
# CL4 — testUnknownExtension: unknown extension (.xyz) exits 1 with extension in message
# ---------------------------------------------------------------------------

timed("CL4", "unknown extension exits 1 with extension in message"):
  let tmpDir = createTempDir("vcfparty_", "")
  let tmpFile = tmpDir / "input.xyz"
  writeFile(tmpFile, "dummy")
  let (outp, code) = run(&"scatter -n 2 -o {tmpDir}/shard {tmpFile}")
  doAssert code != 0, "unknown extension should exit non-zero"
  doAssert ".xyz" in outp, &"error message should contain '.xyz', got: {outp}"
  removeDir(tmpDir)

# ---------------------------------------------------------------------------
# CL5 — testBcfNoIndex: BCF without index exits 1 (no auto-scan fallback)
# ---------------------------------------------------------------------------

timed("CL5", "BCF no index exits 1"):
  doAssert fileExists(SmallBcf), "BCF fixture missing — run generate_fixtures.sh"
  let tmpDir = createTempDir("vcfparty_", "")
  let tmpBcf = tmpDir / "noindex.bcf"
  copyFile(SmallBcf, tmpBcf)
  # Intentionally omit .csi so BCF has no index.
  let (outp, code) = run(&"scatter -n 2 -o {tmpDir}/shard {tmpBcf}")
  doAssert code != 0, &"BCF with no index should exit non-zero, got {code}:\n{outp}"
  removeDir(tmpDir)

# ---------------------------------------------------------------------------
# CL6 — testBcfRunForceScan: --force-scan with BCF via run exits 1
# ---------------------------------------------------------------------------

timed("CL6", "BCF run --force-scan exits 1"):
  doAssert fileExists(SmallBcf), "BCF fixture missing — run generate_fixtures.sh"
  let tmpDir = createTempDir("vcfparty_", "")
  let (outp, code) = run(&"run -n 2 -o {tmpDir}/out.vcf.gz --force-scan {SmallBcf} --- cat")
  doAssert code != 0, "--force-scan with BCF via run should exit non-zero"
  doAssert "force-scan" in outp.toLowerAscii,
    &"error should mention force-scan, got: {outp}"
  removeDir(tmpDir)

# ---------------------------------------------------------------------------
# CL7 — testBcfForceScan: --force-scan with BCF via scatter exits 1
# ---------------------------------------------------------------------------

timed("CL7", "BCF --force-scan exits 1"):
  doAssert fileExists(SmallBcf), "BCF fixture missing — run generate_fixtures.sh"
  let tmpDir = createTempDir("vcfparty_", "")
  let (outp, code) = run(&"scatter -n 2 -o {tmpDir}/shard --force-scan {SmallBcf}")
  doAssert code != 0, "--force-scan with BCF should exit non-zero"
  doAssert "force-scan" in outp.toLowerAscii,
    &"error should mention force-scan, got: {outp}"
  removeDir(tmpDir)

# ===========================================================================
# CL8–CL13 — Integration: scatter correctness for VCF and BCF
# ===========================================================================

# ---------------------------------------------------------------------------
# CL8 — testMissingIndex: no index file → auto-scan fallback with warning; shards produced
# ---------------------------------------------------------------------------

timed("CL8", "no-index scatter: auto-scan with warning, shards produced"):
  # With no index, scatter should warn and fall back to BGZF scan automatically.
  let tmpDir = createTempDir("vcfparty_", "")
  let tmpVcf = tmpDir / "noindex.vcf.gz"
  copyFile(SmallVcf, tmpVcf)
  let outp_template = tmpDir / "out.vcf.gz"
  let (outp, code) = run(&"scatter -n 2 -o {outp_template} {tmpVcf}")
  doAssert code == 0, &"no-index scatter should succeed (auto-scan), got exit {code}:\n{outp}"
  doAssert "warning" in outp.toLowerAscii,
    &"no-index scatter should print a warning, got: {outp}"
  doAssert fileExists(tmpDir / "shard_1.out.vcf.gz"), "shard 1 missing after no-index scatter"
  doAssert fileExists(tmpDir / "shard_2.out.vcf.gz"), "shard 2 missing after no-index scatter"
  removeDir(tmpDir)

# ---------------------------------------------------------------------------
# CL9 — testForceScanFlag: --force-scan ignores existing index; shards valid, record count matches
# ---------------------------------------------------------------------------

timed("CL9", "--force-scan: 4 shards, records match"):
  let tmpDir = createTempDir("vcfparty_", "")
  let outp_template = tmpDir / "out.vcf.gz"
  let (runOutp, runCode) = run(&"scatter -n 4 -o {outp_template} --force-scan {SmallVcf}")
  doAssert runCode == 0, &"--force-scan exited non-zero:\n{runOutp}"
  proc countAndCheckFS(path: string): int =
    ## Count records; also validates the file (bcftools exits 0 on success).
    let (o, code) = execCmdEx("bcftools view -HG " & path & " 2>/dev/null")
    doAssert code == 0, &"bcftools rejected --force-scan shard: {path}"
    o.splitLines.countIt(it.len > 0)
  var total = 0
  for i in 1..4:
    let p = tmpDir / ("shard_" & $i & ".out.vcf.gz")
    doAssert fileExists(p), &"--force-scan shard {i} missing"
    total += countAndCheckFS(p)
  doAssert total == countAndCheckFS(SmallVcf),
    &"--force-scan record count mismatch: shards={total}"
  removeDir(tmpDir)

# ---------------------------------------------------------------------------
# CL10 — testEndToEnd: scatter -n 4 VCF; all shards valid, content hash matches
# ---------------------------------------------------------------------------

timed("CL10", "e2e: record content hash matches original"):
  let tmpDir = createTempDir("vcfparty_", "")
  let outp_template = tmpDir / "out.vcf.gz"

  let (runOutp, runCode) = run(&"scatter -n 4 -o {outp_template} {SmallVcf}")
  doAssert runCode == 0, &"vcfparty scatter exited non-zero:\n{runOutp}"

  for i in 1..4:
    let shardPath = tmpDir / ("shard_" & $i & ".out.vcf.gz")
    doAssert fileExists(shardPath), &"shard {i} not found: {shardPath}"
    let (bcfOutp, bcfCode) = execCmdEx(
      "bcftools view -HG " & shardPath & " > /dev/null 2>&1")
    doAssert bcfCode == 0,
      &"bcftools rejected shard {i}: {bcfOutp}"

  var shardPaths: seq[string]
  for i in 1..4:
    shardPaths.add(tmpDir / ("shard_" & $i & ".out.vcf.gz"))

  doAssert recordsHash(shardPaths) == recordsHash(@[SmallVcf]),
    "e2e content hash mismatch: record corruption or reordering detected"

  removeDir(tmpDir)

# ---------------------------------------------------------------------------
# CL11 — testCsiIndex: CSI-indexed VCF scatter; shards valid, record count matches
# ---------------------------------------------------------------------------

timed("CL11", "CSI: scatter -n 4, all shards valid"):
  doAssert fileExists(CsiVcf), "CSI fixture missing — run generate_fixtures.sh"
  let tmpDir = createTempDir("vcfparty_", "")
  let outp_template = tmpDir / "out.vcf.gz"

  let (runOutp, runCode) = run(&"scatter -n 4 -o {outp_template} {CsiVcf}")
  doAssert runCode == 0, &"vcfparty scatter (CSI) exited non-zero:\n{runOutp}"

  proc countRec(path: string): int =
    let (o, _) = execCmdEx("bcftools view -HG " & path & " 2>/dev/null | wc -l")
    o.strip.parseInt

  var shardTotal = 0
  for i in 1..4:
    let path = tmpDir / ("shard_" & $i & ".out.vcf.gz")
    doAssert fileExists(path), &"CSI shard {i} not found: {path}"
    let (bcfO, bcfC) = execCmdEx("bcftools view -HG " & path & " > /dev/null 2>&1")
    doAssert bcfC == 0, &"bcftools rejected CSI shard {i}: {bcfO}"
    shardTotal += countRec(path)
  let origTotal = countRec(CsiVcf)
  doAssert shardTotal == origTotal,
    &"CSI record count mismatch: shards={shardTotal} orig={origTotal}"

  removeDir(tmpDir)

# ---------------------------------------------------------------------------
# CL12 — testBcfExtension: BCF scatter produces .bcf shards with matching content hash
# ---------------------------------------------------------------------------

timed("CL12", "BCF .bcf extension -> BCF mode, content hash matches"):
  doAssert fileExists(SmallBcf), "BCF fixture missing — run generate_fixtures.sh"
  let tmpDir = createTempDir("vcfparty_", "")
  let outp_template = tmpDir / "out.bcf"
  let (outp, code) = run(&"scatter -n 4 -o {outp_template} {SmallBcf}")
  doAssert code == 0, &"BCF scatter exited non-zero:\n{outp}"
  var bcfShardPaths: seq[string]
  for i in 1..4:
    let p = tmpDir / ("shard_" & $i & ".out.bcf")
    doAssert fileExists(p), &"BCF shard {i} missing (.bcf extension): {p}"
    let (bo, bc) = execCmdEx("bcftools view -HG " & p & " > /dev/null 2>&1")
    doAssert bc == 0, &"bcftools rejected BCF shard {i}: {bo}"
    bcfShardPaths.add(p)
  doAssert recordsHash(bcfShardPaths) == recordsHash(@[SmallBcf]),
    "BCF content hash mismatch: record corruption or reordering detected"
  removeDir(tmpDir)

# ---------------------------------------------------------------------------
# CL13 — testKg1000Genomes: large 1KG VCF scatter -n 10 (skipped if fixture absent)
# ---------------------------------------------------------------------------

block testKg1000Genomes:
  if not fileExists(KgVcf):
    echo "SKIP CL13 1KG chr22: file not present (run tests/generate_fixtures.sh)"
  else:
    timed("CL13", "1KG chr22 scatter: records match original"):
      let tmpDir = createTempDir("vcfparty_", "")
      let outp_template = tmpDir / "out.vcf.gz"

      let (runOutp, runCode) = run(&"scatter -n 10 -o {outp_template} {KgVcf}")
      doAssert runCode == 0, &"vcfparty scatter (1KG) exited non-zero:\n{runOutp}"

      # With -n 10, nDigits=2 so names are shard_01.out.vcf.gz … shard_10.out.vcf.gz
      for i in 1..10:
        let path = tmpDir / ("shard_" & align($i, 2, '0') & ".out.vcf.gz")
        doAssert fileExists(path), &"1KG shard {i} not found: {path}"
        let (bcfOutp, bcfCode) = execCmdEx(
          "bcftools view -HG " & path & " > /dev/null 2>&1")
        doAssert bcfCode == 0, &"bcftools rejected 1KG shard {i}: {bcfOutp}"

      proc countRecs(path: string): int =
        let (o, _) = execCmdEx("bcftools view -HG " & path & " 2>/dev/null | wc -l")
        o.strip.parseInt

      var shardTotal = 0
      for i in 1..10:
        shardTotal += countRecs(tmpDir / ("shard_" & align($i, 2, '0') & ".out.vcf.gz"))
      let origTotal = countRecs(KgVcf)
      doAssert shardTotal == origTotal,
        &"1KG record count mismatch: shards={shardTotal} orig={origTotal}"

      removeDir(tmpDir)

# ===========================================================================
# CL30 — -u flag (force uncompressed output)
# ===========================================================================

proc isBgzf(path: string): bool =
  let f = open(path, fmRead)
  var b: array[2, byte]
  discard f.readBytes(b, 0, 2)
  f.close()
  b[0] == 0x1f'u8 and b[1] == 0x8b'u8

# ---------------------------------------------------------------------------
# CL30.1 — run -u +concat+: output is uncompressed despite .vcf.gz extension
# ---------------------------------------------------------------------------

timed("CL30.1", "run -u +concat+: uncompressed output, warning emitted"):
  let tmpDir = createTempDir("vcfparty_", "")
  let outFile = tmpDir / "out.vcf.gz"
  let (outp, code) = run(
    &"run -n 2 -u -o {outFile} {SmallVcf} ::: cat +concat+")
  doAssert code == 0, &"CL30.1 run -u +concat+ exited {code}:\n{outp}"
  doAssert fileExists(outFile), "CL30.1: output missing"
  doAssert not isBgzf(outFile), "CL30.1: output should be uncompressed despite .vcf.gz extension"
  doAssert "warning" in outp.toLowerAscii,
    &"CL30.1: expected extension mismatch warning, got:\n{outp}"
  removeDir(tmpDir)

# ---------------------------------------------------------------------------
# CL30.2 — run -u with {}: exits non-zero
# ---------------------------------------------------------------------------

timed("CL30.2", "run -u with {}: exits non-zero"):
  let (outp, code) = run(
    &"run -n 2 -u {SmallVcf} ::: cat -o out.{{}}.vcf.gz")
  doAssert code != 0, &"CL30.2: -u with {{}} should exit non-zero, got {code}"
  doAssert "tool-managed" in outp.toLowerAscii or "{}" in outp,
    &"CL30.2: expected tool-managed error, got:\n{outp}"

# ---------------------------------------------------------------------------
# CL30.3 — run -u +concat+ with .vcf extension: no warning, uncompressed
# ---------------------------------------------------------------------------

timed("CL30.3", "run -u +concat+ .vcf: uncompressed, no warning"):
  let tmpDir = createTempDir("vcfparty_", "")
  let outFile = tmpDir / "out.vcf"
  let (outp, code) = run(
    &"run -n 2 -u -o {outFile} {SmallVcf} ::: cat +concat+")
  doAssert code == 0, &"CL30.3 run -u +concat+ .vcf exited {code}:\n{outp}"
  doAssert fileExists(outFile), "CL30.3: output missing"
  doAssert not isBgzf(outFile), "CL30.3: output should be uncompressed"
  let cnt = execCmdEx("grep -c '^[^#]' " & outFile & " 2>/dev/null || true")[0].strip.parseInt
  doAssert cnt > 0, "CL30.3: output has no records"
  removeDir(tmpDir)

# ===========================================================================
# CL27 — +merge+ basic integration
# ===========================================================================

timed("CL27", "+merge+: exits zero, records present"):
  let tmpDir  = createTempDir("vcfparty_", "")
  let outFile = tmpDir / "out.vcf"
  let (outp, code) = run(
    &"run -n 2 -o {outFile} {SmallVcf} ::: bcftools view -Ov +merge+")
  doAssert code == 0, &"C27: +merge+ exited {code}:\n{outp}"
  doAssert fileExists(outFile), "C27: output file missing"
  let origCnt = execCmdEx("bcftools view -H " & SmallVcf & " 2>/dev/null | wc -l")[0].strip.parseInt
  var outCnt = 0
  for line in lines(outFile):
    if not line.startsWith("#"): inc outCnt
  doAssert outCnt == origCnt,
    &"C27: record count mismatch: orig={origCnt} out={outCnt}"
  removeDir(tmpDir)

# ===========================================================================
# CL28 — FD inheritance regression: +concat+ with slow subprocess
# ===========================================================================
#
# Regression test for the fd-inheritance deadlock fixed by FD_CLOEXEC.
#
# Without FD_CLOEXEC on stdinPipe[1]: when child i+1 is forked it inherits
# stdinPipe_i[1] from the writer thread. After the writer closes its copy,
# child i's subprocess (cat reading stdin) never sees EOF because child i+1
# still holds the fd. Child i hangs → interceptor i hangs → +concat+ deadlocks.
#
# The sleep before cat ensures the writer for shard i is still alive when
# shard i+1 is forked, maximising the chance of the race being triggered.

timed("CL28", "+concat+ fd-inheritance: 4 shards, slow subprocess"):
  let tmpDir  = createTempDir("vcfparty_", "")
  let outFile = tmpDir / "out.vcf.gz"
  # 4 shards, slow subprocess: sleep then cat (passes data through unchanged).
  let (outp, code) = run(
    &"run -n 4 -o {outFile} {SmallVcf} ::: sh -c 'sleep 0.05 && cat' +concat+")
  doAssert code == 0,
    &"C28: +concat+ with slow subprocess exited {code}:\n{outp}"
  doAssert fileExists(outFile), "C28: output file missing"
  let origCnt = execCmdEx(
    "bcftools view -H " & SmallVcf & " 2>/dev/null | wc -l")[0].strip.parseInt
  let outCnt  = execCmdEx(
    "bcftools view -H " & outFile  & " 2>/dev/null | wc -l")[0].strip.parseInt
  doAssert outCnt == origCnt,
    &"C28: record count mismatch: orig={origCnt} out={outCnt}"
  removeDir(tmpDir)

# ===========================================================================
# CL31 — -n exceeds index entry count
# ===========================================================================
# Scatter and run paths must reject nShards > available index entries with a
# clear error message rather than producing empty shards. Tests cover:
#   CL31.1 — scatter VCF: error suggests --force-scan
#   CL31.2 — scatter BCF: error does NOT suggest --force-scan (BCF can't scan)
#   CL31.3 — scatter VCF --force-scan: succeeds when raw blocks > requested n
#   CL31.4 — run +merge+ VCF: same error path as scatter
#   CL31.5 — -n at the exact limit succeeds (boundary check)

# ---------------------------------------------------------------------------
# CL31.1 — scatter -n way too high on VCF: errors, suggests --force-scan
# ---------------------------------------------------------------------------

timed("CL31.1", "scatter -n too high VCF: errors with --force-scan suggestion"):
  let tmpDir = createTempDir("vcfparty_", "")
  # chr22_1kg.vcf.gz has only 57 index entries; -n 100 must error.
  let (outp, code) = run(&"scatter -n 100 -o {tmpDir}/shard_{{}}.vcf.gz {KgVcf}")
  doAssert code != 0,
    &"CL31.1: -n 100 on chr22_1kg.vcf.gz (57 voffs) should error, got {code}:\n{outp}"
  doAssert "index entries" in outp,
    &"CL31.1: error should mention 'index entries', got: {outp}"
  doAssert "--force-scan" in outp,
    &"CL31.1: VCF error should suggest --force-scan, got: {outp}"
  removeDir(tmpDir)

# ---------------------------------------------------------------------------
# CL31.2 — scatter -n way too high on BCF: errors, NO --force-scan suggestion
# ---------------------------------------------------------------------------

timed("CL31.2", "scatter -n too high BCF: errors without --force-scan suggestion"):
  doAssert fileExists(DataDir / "chr22_1kg.bcf"),
    "chr22_1kg.bcf fixture missing — run generate_fixtures.sh"
  let tmpDir = createTempDir("vcfparty_", "")
  # chr22_1kg.bcf has only 58 index entries; -n 100 must error.
  let (outp, code) = run(
    &"scatter -n 100 -o {tmpDir}/shard_{{}}.bcf {DataDir}/chr22_1kg.bcf")
  doAssert code != 0,
    &"CL31.2: -n 100 on chr22_1kg.bcf (58 voffs) should error, got {code}:\n{outp}"
  doAssert "index entries" in outp,
    &"CL31.2: error should mention 'index entries', got: {outp}"
  doAssert "--force-scan" notin outp,
    &"CL31.2: BCF error should NOT suggest --force-scan, got: {outp}"
  removeDir(tmpDir)

# ---------------------------------------------------------------------------
# CL31.3 — scatter --force-scan VCF: succeeds when raw blocks exceed -n
# ---------------------------------------------------------------------------

timed("CL31.3", "scatter --force-scan recovers: 100 shards, records match"):
  let tmpDir = createTempDir("vcfparty_", "")
  # chr22_1kg.vcf.gz has 57 voffs but ~834 raw BGZF blocks; --force-scan
  # gives access to all raw blocks so -n 100 works.
  let (outp, code) = run(
    &"scatter -n 100 --force-scan -o {tmpDir}/shard_{{}}.vcf.gz {KgVcf}")
  doAssert code == 0,
    &"CL31.3: --force-scan with -n 100 should succeed, got {code}:\n{outp}"
  let shards = toSeq(walkFiles(tmpDir / "shard_*.vcf.gz"))
  doAssert shards.len == 100,
    &"CL31.3: expected 100 shards, got {shards.len}"
  # Sanity-check record count
  let origCnt = execCmdEx(
    "bcftools view -H " & KgVcf & " 2>/dev/null | wc -l")[0].strip.parseInt
  var outCnt = 0
  for s in shards:
    outCnt += execCmdEx(
      "bcftools view -H " & s & " 2>/dev/null | wc -l")[0].strip.parseInt
  doAssert outCnt == origCnt,
    &"CL31.3: record count mismatch: orig={origCnt} shards={outCnt}"
  removeDir(tmpDir)

# ---------------------------------------------------------------------------
# CL31.4 — run +merge+ -n too high errors with same path
# ---------------------------------------------------------------------------

timed("CL31.4", "run +merge+ -n too high: errors"):
  doAssert fileExists(DataDir / "chr22_1kg.bcf"),
    "chr22_1kg.bcf fixture missing — run generate_fixtures.sh"
  let tmpDir = createTempDir("vcfparty_", "")
  let (outp, code) = run(
    &"run -n 100 -o {tmpDir}/out.vcf.gz {DataDir}/chr22_1kg.bcf ::: bcftools view -Ov +merge+")
  doAssert code != 0,
    &"CL31.4: run +merge+ -n 100 on chr22_1kg.bcf should error, got {code}:\n{outp}"
  doAssert "index entries" in outp,
    &"CL31.4: error should mention 'index entries', got: {outp}"
  removeDir(tmpDir)

# ---------------------------------------------------------------------------
# CL31.5 — -n at the exact voff count succeeds
# ---------------------------------------------------------------------------

timed("CL31.5", "scatter -n at exact voff limit: 28 shards"):
  let tmpDir = createTempDir("vcfparty_", "")
  # chr22_1kg.vcf.gz has 28 voffs; -n 28 must succeed.
  let (outp, code) = run(&"scatter -n 28 -o {tmpDir}/shard_{{}}.vcf.gz {KgVcf}")
  doAssert code == 0,
    &"CL31.5: -n 28 (== voff count) should succeed, got {code}:\n{outp}"
  let shards = toSeq(walkFiles(tmpDir / "shard_*.vcf.gz"))
  doAssert shards.len == 28,
    &"CL31.5: expected 28 shards, got {shards.len}"
  removeDir(tmpDir)

# ===========================================================================
# CL32 — --clamp-shards reduces -n instead of erroring
# ===========================================================================
# When --clamp-shards is set and -n exceeds the available index entries, the
# tool prints an info line and silently produces fewer shards rather than
# erroring out (CL31's no-clamp behaviour).

# ---------------------------------------------------------------------------
# CL32.1 — scatter --clamp-shards on VCF: -n 100 → 57 shards
# ---------------------------------------------------------------------------

timed("CL32.1", "scatter --clamp-shards VCF: 28 shards, records match"):
  let tmpDir = createTempDir("vcfparty_", "")
  # chr22_1kg.vcf.gz has 28 index entries; -n 100 must clamp to 28.
  let (outp, code) = run(
    &"scatter -n 100 --clamp-shards -o {tmpDir}/shard_{{}}.vcf.gz {KgVcf}")
  doAssert code == 0,
    &"CL32.1: --clamp-shards should succeed, got {code}:\n{outp}"
  doAssert "clamp-shards" in outp,
    &"CL32.1: info should mention 'clamp-shards', got: {outp}"
  let shards = toSeq(walkFiles(tmpDir / "shard_*.vcf.gz"))
  doAssert shards.len == 28,
    &"CL32.1: expected 28 clamped shards, got {shards.len}"
  # Sanity-check record count matches the original.
  let origCnt = execCmdEx(
    "bcftools view -H " & KgVcf & " 2>/dev/null | wc -l")[0].strip.parseInt
  var outCnt = 0
  for s in shards:
    outCnt += execCmdEx(
      "bcftools view -H " & s & " 2>/dev/null | wc -l")[0].strip.parseInt
  doAssert outCnt == origCnt,
    &"CL32.1: record count mismatch: orig={origCnt} shards={outCnt}"
  removeDir(tmpDir)


# ===========================================================================
# CL33 — streaming +concat+ correctness on a non-trivial fixture
# ===========================================================================
# After the runInterceptor refactor from buffer-then-process to streaming, the
# +concat+ output must remain byte-/record-equivalent across all (input,
# output) compression combinations on a real ~25k-record fixture, and must
# fail cleanly on #CHROM mismatch and on header-only shards.

# ---------------------------------------------------------------------------
# CL33.1 — +concat+ BGZF input → BGZF output (the dominant production path)
# ---------------------------------------------------------------------------

timed("CL33.1", "streaming +concat+ BGZF->BGZF: records, hash match"):
  let tmpDir = createTempDir("vcfparty_", "")
  let outPath = tmpDir / "out.vcf.gz"
  let (outp, code) = run(
    &"run -n 4 -o {outPath} {KgVcf} ::: bcftools view -Oz +concat+")
  doAssert code == 0,
    &"CL33.1: +concat+ BGZF→BGZF on chr22_1kg should succeed, got {code}:\n{outp}"
  let origCnt = execCmdEx(
    "bcftools view -H " & KgVcf & " 2>/dev/null | wc -l")[0].strip.parseInt
  let outCnt = execCmdEx(
    "bcftools view -H " & outPath & " 2>/dev/null | wc -l")[0].strip.parseInt
  doAssert outCnt == origCnt,
    &"CL33.1: record count mismatch: orig={origCnt} concat={outCnt}"
  let origHash = recordsHash(@[KgVcf])
  let outHash  = recordsHash(@[outPath])
  doAssert origHash == outHash,
    &"CL33.1: record content hash mismatch:\n  orig={origHash}\n  out ={outHash}"
  removeDir(tmpDir)

# ---------------------------------------------------------------------------
# CL33.2 — +concat+ BGZF input → uncompressed output (-u)
# ---------------------------------------------------------------------------

timed("CL33.2", "streaming +concat+ BGZF->uncompressed: records match"):
  let tmpDir = createTempDir("vcfparty_", "")
  let outPath = tmpDir / "out.vcf"
  let (outp, code) = run(
    &"run -n 4 -u -o {outPath} {KgVcf} ::: bcftools view -Oz +concat+")
  doAssert code == 0,
    &"CL33.2: +concat+ -u should succeed, got {code}:\n{outp}"
  let origCnt = execCmdEx(
    "bcftools view -H " & KgVcf & " 2>/dev/null | wc -l")[0].strip.parseInt
  let outCnt = execCmdEx(
    "bcftools view -H " & outPath & " 2>/dev/null | wc -l")[0].strip.parseInt
  doAssert outCnt == origCnt,
    &"CL33.2: record count mismatch: orig={origCnt} concat={outCnt}"
  let origHash = recordsHash(@[KgVcf])
  let outHash  = recordsHash(@[outPath])
  doAssert origHash == outHash,
    &"CL33.2: record content hash mismatch"
  removeDir(tmpDir)

# ---------------------------------------------------------------------------
# CL33.3 — +concat+ uncompressed pipeline → BGZF output (recompression mode)
# ---------------------------------------------------------------------------

timed("CL33.3", "streaming +concat+ uncompressed->BGZF: records match"):
  let tmpDir = createTempDir("vcfparty_", "")
  let outPath = tmpDir / "out.vcf.gz"
  # bcftools view -Ov emits uncompressed VCF; vcfparty must recompress to BGZF.
  let (outp, code) = run(
    &"run -n 4 -o {outPath} {KgVcf} ::: bcftools view -Ov +concat+")
  doAssert code == 0,
    &"CL33.3: +concat+ uncompressed→BGZF should succeed, got {code}:\n{outp}"
  let origCnt = execCmdEx(
    "bcftools view -H " & KgVcf & " 2>/dev/null | wc -l")[0].strip.parseInt
  let outCnt = execCmdEx(
    "bcftools view -H " & outPath & " 2>/dev/null | wc -l")[0].strip.parseInt
  doAssert outCnt == origCnt,
    &"CL33.3: record count mismatch: orig={origCnt} concat={outCnt}"
  let origHash = recordsHash(@[KgVcf])
  let outHash  = recordsHash(@[outPath])
  doAssert origHash == outHash,
    &"CL33.3: record content hash mismatch"
  removeDir(tmpDir)

# ---------------------------------------------------------------------------
# CL33.4 — #CHROM mismatch on the streaming +concat+ path fails clean
# ---------------------------------------------------------------------------
# Use sed to rewrite a sample column in the bcftools output of one shard
# (different shards via {} substitution).  Shard 1 keeps the original #CHROM,
# shard 2 has a renamed sample column, and the streaming validation must
# detect the mismatch and exit non-zero.

timed("CL33.4", "streaming +concat+ #CHROM mismatch: clean failure"):
  let tmpDir = createTempDir("vcfparty_", "")
  # Rename one sample column (HG00096 → ZZ00000) in only shard 2; all other
  # shards pass through unchanged.  bcftools accepts the rewritten header
  # because it is still a valid #CHROM line — the only thing that differs is
  # one sample name — which is exactly what vcfparty's #CHROM validation must
  # detect.  bcftools reheader is the canonical tool for this.
  let renameFile = tmpDir / "rename.txt"
  writeFile(renameFile, "HG00096 ZZ00000\n")
  let cmd = &"run -n 4 -o {tmpDir}/out.vcf.gz {KgVcf} ::: " &
            &"sh -c 'if [ {{}} = 2 ]; then bcftools reheader -s {renameFile} | bcftools view -Oz; else bcftools view -Oz; fi' " &
            "+concat+"
  let (outp, code) = run(cmd)
  doAssert code != 0,
    &"CL33.4: #CHROM mismatch should fail, got {code}:\n{outp}"
  doAssert "CHROM" in outp,
    &"CL33.4: error should mention #CHROM, got: {outp}"
  removeDir(tmpDir)

# ---------------------------------------------------------------------------
# CL33.5 — header-only shard handled correctly
# ---------------------------------------------------------------------------
# Subprocess `bcftools view -h` emits only the header (no records) for shard 1
# only; other shards emit full records.  Streaming must handle the case where
# phase B's header detection consumes the entire stream and phase C runs on an
# empty record buffer.

timed("CL33.5", "streaming +concat+ header-only shard: handled correctly"):
  let tmpDir = createTempDir("vcfparty_", "")
  let outPath = tmpDir / "out.vcf.gz"
  let cmd = &"run -n 4 -o {outPath} {KgVcf} ::: " &
            "sh -c 'if [ {} = 1 ]; then bcftools view -h -Oz; else bcftools view -Oz; fi' " &
            "+concat+"
  let (outp, code) = run(cmd)
  doAssert code == 0,
    &"CL33.5: header-only shard should not crash, got {code}:\n{outp}"
  let outCnt = execCmdEx(
    "bcftools view -H " & outPath & " 2>/dev/null | wc -l")[0].strip.parseInt
  # Shard 1 contributed zero records; shards 2..4 contributed their full
  # share.  We don't assert the exact remaining count (it depends on shard
  # boundaries) but assert it is non-zero and strictly less than the original.
  let origCnt = execCmdEx(
    "bcftools view -H " & KgVcf & " 2>/dev/null | wc -l")[0].strip.parseInt
  doAssert outCnt > 0 and outCnt < origCnt,
    &"CL33.5: expected 0 < outCnt < origCnt, got outCnt={outCnt} origCnt={origCnt}"
  removeDir(tmpDir)
