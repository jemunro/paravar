## Tests for CLI argument handling.
## Run from project root: nim c -r tests/test_cli.nim
## Requires the vcfparty binary to be built first (nimble build).

import std/[os, osproc, sequtils, strformat, strutils, tempfiles]

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

block buildBinary:
  if not fileExists(BinPath):
    let (outp, code) = execCmdEx("nimble build 2>&1")
    if code != 0:
      echo "nimble build failed:\n", outp
      quit(1)
  doAssert fileExists(BinPath), "binary not found: " & BinPath & " (run nimble build)"
  echo "PASS CL0 binary available"

# ===========================================================================
# CL1–CL7 — Error cases (missing flags, invalid args, unsupported modes)
# ===========================================================================

# ---------------------------------------------------------------------------
# CL1 — testMissingN: missing -n exits non-zero with 'n' in error message
# ---------------------------------------------------------------------------

block testMissingN:
  let (outp, code) = run(&"scatter -o /tmp/vcfparty_cli_test {SmallVcf}")
  doAssert code != 0, "missing -n should exit non-zero"
  doAssert "n" in outp.toLowerAscii,
    &"missing -n error should mention 'n', got: {outp}"
  echo "PASS CL1 missing -n exits non-zero"

# ---------------------------------------------------------------------------
# CL2 — testInvalidN0: -n 0 exits non-zero
# ---------------------------------------------------------------------------

block testInvalidN0:
  let (_, code) = run(&"scatter -n 0 -o /tmp/vcfparty_cli_test {SmallVcf}")
  doAssert code != 0, "-n 0 should exit non-zero"
  echo "PASS CL2 -n 0 exits non-zero"

# ---------------------------------------------------------------------------
# CL3 — testMissingO: missing -o exits non-zero with 'o' in error message
# ---------------------------------------------------------------------------

block testMissingO:
  let (outp, code) = run(&"scatter -n 2 {SmallVcf}")
  doAssert code != 0, "missing -o should exit non-zero"
  doAssert "o" in outp.toLowerAscii,
    &"missing -o error should mention 'o', got: {outp}"
  echo "PASS CL3 missing -o exits non-zero"

# ---------------------------------------------------------------------------
# CL4 — testUnknownExtension: unknown extension (.xyz) exits 1 with extension in message
# ---------------------------------------------------------------------------

block testUnknownExtension:
  let tmpDir = createTempDir("vcfparty_", "")
  let tmpFile = tmpDir / "input.xyz"
  writeFile(tmpFile, "dummy")
  let (outp, code) = run(&"scatter -n 2 -o {tmpDir}/shard {tmpFile}")
  doAssert code != 0, "unknown extension should exit non-zero"
  doAssert ".xyz" in outp, &"error message should contain '.xyz', got: {outp}"
  removeDir(tmpDir)
  echo "PASS CL4 unknown extension exits 1 with extension in message"

# ---------------------------------------------------------------------------
# CL5 — testBcfNoIndex: BCF without index exits 1 (no auto-scan fallback)
# ---------------------------------------------------------------------------

block testBcfNoIndex:
  doAssert fileExists(SmallBcf), "BCF fixture missing — run generate_fixtures.sh"
  let tmpDir = createTempDir("vcfparty_", "")
  let tmpBcf = tmpDir / "noindex.bcf"
  copyFile(SmallBcf, tmpBcf)
  # Intentionally omit .csi so BCF has no index.
  let (outp, code) = run(&"scatter -n 2 -o {tmpDir}/shard {tmpBcf}")
  doAssert code != 0, &"BCF with no index should exit non-zero, got {code}:\n{outp}"
  removeDir(tmpDir)
  echo "PASS CL5 BCF no index exits 1"

# ---------------------------------------------------------------------------
# CL6 — testBcfRunForceScan: --force-scan with BCF via run exits 1
# ---------------------------------------------------------------------------

block testBcfRunForceScan:
  doAssert fileExists(SmallBcf), "BCF fixture missing — run generate_fixtures.sh"
  let tmpDir = createTempDir("vcfparty_", "")
  let (outp, code) = run(&"run -n 2 -o {tmpDir}/out.vcf.gz --force-scan {SmallBcf} --- cat")
  doAssert code != 0, "--force-scan with BCF via run should exit non-zero"
  doAssert "force-scan" in outp.toLowerAscii,
    &"error should mention force-scan, got: {outp}"
  removeDir(tmpDir)
  echo "PASS CL6 BCF run --force-scan exits 1"

# ---------------------------------------------------------------------------
# CL7 — testBcfForceScan: --force-scan with BCF via scatter exits 1
# ---------------------------------------------------------------------------

block testBcfForceScan:
  doAssert fileExists(SmallBcf), "BCF fixture missing — run generate_fixtures.sh"
  let tmpDir = createTempDir("vcfparty_", "")
  let (outp, code) = run(&"scatter -n 2 -o {tmpDir}/shard --force-scan {SmallBcf}")
  doAssert code != 0, "--force-scan with BCF should exit non-zero"
  doAssert "force-scan" in outp.toLowerAscii,
    &"error should mention force-scan, got: {outp}"
  removeDir(tmpDir)
  echo "PASS CL7 BCF --force-scan exits 1"

# ===========================================================================
# CL8–CL13 — Integration: scatter correctness for VCF and BCF
# ===========================================================================

# ---------------------------------------------------------------------------
# CL8 — testMissingIndex: no index file → auto-scan fallback with warning; shards produced
# ---------------------------------------------------------------------------

block testMissingIndex:
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
  echo "PASS CL8 no-index scatter: auto-scan with warning, shards produced"

# ---------------------------------------------------------------------------
# CL9 — testForceScanFlag: --force-scan ignores existing index; shards valid, record count matches
# ---------------------------------------------------------------------------

block testForceScanFlag:
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
  echo &"PASS CL9 --force-scan: 4 shards, {total} records"
  removeDir(tmpDir)

# ---------------------------------------------------------------------------
# CL10 — testEndToEnd: scatter -n 4 VCF; all shards valid, content hash matches
# ---------------------------------------------------------------------------

block testEndToEnd:
  let tmpDir = createTempDir("vcfparty_", "")
  let outp_template = tmpDir / "out.vcf.gz"

  let (runOutp, runCode) = run(&"scatter -n 4 -o {outp_template} {SmallVcf}")
  doAssert runCode == 0, &"vcfparty scatter exited non-zero:\n{runOutp}"
  echo "PASS CL10.1 e2e: vcfparty scatter -n 4 exited 0"

  for i in 1..4:
    let shardPath = tmpDir / ("shard_" & $i & ".out.vcf.gz")
    doAssert fileExists(shardPath), &"shard {i} not found: {shardPath}"
    let (bcfOutp, bcfCode) = execCmdEx(
      "bcftools view -HG " & shardPath & " > /dev/null 2>&1")
    doAssert bcfCode == 0,
      &"bcftools rejected shard {i}: {bcfOutp}"
  echo "PASS CL10.2 e2e: all 4 shards are valid VCFs (bcftools view)"

  var shardPaths: seq[string]
  for i in 1..4:
    shardPaths.add(tmpDir / ("shard_" & $i & ".out.vcf.gz"))

  doAssert recordsHash(shardPaths) == recordsHash(@[SmallVcf]),
    "e2e content hash mismatch: record corruption or reordering detected"
  echo "PASS CL10.3 e2e: record content hash matches original (no corruption, correct count)"

  removeDir(tmpDir)

# ---------------------------------------------------------------------------
# CL11 — testCsiIndex: CSI-indexed VCF scatter; shards valid, record count matches
# ---------------------------------------------------------------------------

block testCsiIndex:
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
  echo &"PASS CL11 CSI: scatter -n 4, all shards valid, {origTotal} records"

  removeDir(tmpDir)

# ---------------------------------------------------------------------------
# CL12 — testBcfExtension: BCF scatter produces .bcf shards with matching content hash
# ---------------------------------------------------------------------------

block testBcfExtension:
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
  echo "PASS CL12 BCF .bcf extension → BCF mode, shards have .bcf extension, content hash matches"

# ---------------------------------------------------------------------------
# CL13 — testKg1000Genomes: large 1KG VCF scatter -n 10 (skipped if fixture absent)
# ---------------------------------------------------------------------------

block testKg1000Genomes:
  if not fileExists(KgVcf):
    echo "SKIP CL13 1KG chr22: file not present (run tests/generate_fixtures.sh)"
  else:
    let tmpDir = createTempDir("vcfparty_", "")
    let outp_template = tmpDir / "out.vcf.gz"

    let (runOutp, runCode) = run(&"scatter -n 10 -o {outp_template} {KgVcf}")
    doAssert runCode == 0, &"vcfparty scatter (1KG) exited non-zero:\n{runOutp}"
    echo "PASS CL13.1 1KG: vcfparty scatter -n 10 exited 0"

    # With -n 10, nDigits=2 so names are shard_01.out.vcf.gz … shard_10.out.vcf.gz
    for i in 1..10:
      let path = tmpDir / ("shard_" & align($i, 2, '0') & ".out.vcf.gz")
      doAssert fileExists(path), &"1KG shard {i} not found: {path}"
      let (bcfOutp, bcfCode) = execCmdEx(
        "bcftools view -HG " & path & " > /dev/null 2>&1")
      doAssert bcfCode == 0, &"bcftools rejected 1KG shard {i}: {bcfOutp}"
    echo "PASS CL13.2 1KG: all 10 shards valid VCFs (bcftools view)"

    proc countRecs(path: string): int =
      let (o, _) = execCmdEx("bcftools view -HG " & path & " 2>/dev/null | wc -l")
      o.strip.parseInt

    var shardTotal = 0
    for i in 1..10:
      shardTotal += countRecs(tmpDir / ("shard_" & align($i, 2, '0') & ".out.vcf.gz"))
    let origTotal = countRecs(KgVcf)
    doAssert shardTotal == origTotal,
      &"1KG record count mismatch: shards={shardTotal} orig={origTotal}"
    echo &"PASS CL13.3 1KG: record count matches original ({origTotal} records across 10 shards)"

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

block testRunUncompressConcat:
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
  echo "PASS CL30.1 run -u +concat+: uncompressed output, warning emitted"

# ---------------------------------------------------------------------------
# CL30.2 — run -u with {}: exits non-zero
# ---------------------------------------------------------------------------

block testRunUncompressBraceError:
  let (outp, code) = run(
    &"run -n 2 -u {SmallVcf} ::: cat -o out.{{}}.vcf.gz")
  doAssert code != 0, &"CL30.2: -u with {{}} should exit non-zero, got {code}"
  doAssert "tool-managed" in outp.toLowerAscii or "{}" in outp,
    &"CL30.2: expected tool-managed error, got:\n{outp}"
  echo "PASS CL30.2 run -u with {}: exits non-zero"

# ---------------------------------------------------------------------------
# CL30.3 — run -u +concat+ with .vcf extension: no warning, uncompressed
# ---------------------------------------------------------------------------

block testRunUncompressConcatVcfExt:
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
  echo "PASS CL30.3 run -u +concat+ .vcf: uncompressed, no warning"

# ===========================================================================
# CL27 — +merge+ basic integration
# ===========================================================================

block testMergeBasic:
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
  echo &"PASS CL27 +merge+: exits zero, {outCnt} records present"

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

block testConcatFdInheritance:
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
  echo &"PASS CL28 +concat+ fd-inheritance: 4 shards, slow subprocess, {outCnt} records"



