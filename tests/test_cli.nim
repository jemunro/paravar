## Tests for CLI argument handling.
## Run from project root: nim c -r tests/test_cli.nim
## Requires the blocky binary to be built first (nimble build).

echo "--------------- Test CLI ---------------"

import std/[os, osproc, sequtils, strformat, strutils, tempfiles]
import test_utils

const BinPath  = "./blocky"
const DataDir  = "tests/data"
const SmallVcf = DataDir / "small.vcf.gz"       # TBI indexed
const CsiVcf   = DataDir / "small_csi.vcf.gz"   # CSI indexed only
const KgVcf    = DataDir / "chr22_1kg.vcf.gz"
const SmallBcf = DataDir / "small.bcf"          # CSI indexed BCF

proc run(args: string): (string, int) =
  ## Run blocky with shell args; combine stdout+stderr; return (outp, code).
  let t = getEnv("BLOCKY_TEST_TIMEOUT", "10")
  execCmdEx("timeout " & t & " " & BinPath & " " & args & " 2>&1")

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
  let (outp, code) = run(&"scatter -o /tmp/blocky_cli_test {SmallVcf}")
  doAssert code != 0, "missing -n should exit non-zero"
  doAssert "n" in outp.toLowerAscii,
    &"missing -n error should mention 'n', got: {outp}"

# ---------------------------------------------------------------------------
# CL2 — testInvalidN0: -n 0 exits non-zero
# ---------------------------------------------------------------------------

timed("CL2", "-n 0 exits non-zero"):
  let (_, code) = run(&"scatter -n 0 -o /tmp/blocky_cli_test {SmallVcf}")
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
  let tmpDir = createTempDir("blocky_", "")
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
  let tmpDir = createTempDir("blocky_", "")
  let tmpBcf = tmpDir / "noindex.bcf"
  copyFile(SmallBcf, tmpBcf)
  # Intentionally omit .csi so BCF has no index.
  let (outp, code) = run(&"scatter -n 2 -o {tmpDir}/shard {tmpBcf}")
  doAssert code != 0, &"BCF with no index should exit non-zero, got {code}:\n{outp}"
  removeDir(tmpDir)

# ---------------------------------------------------------------------------
# CL6 — testBcfRunScan: --scan with BCF via run exits 1
# ---------------------------------------------------------------------------

timed("CL6", "BCF run --scan exits 1"):
  doAssert fileExists(SmallBcf), "BCF fixture missing — run generate_fixtures.sh"
  let tmpDir = createTempDir("blocky_", "")
  let (outp, code) = run(&"run -n 2 -o {tmpDir}/out.vcf.gz --scan {SmallBcf} --- cat")
  doAssert code != 0, "--scan with BCF via run should exit non-zero"
  doAssert "scan" in outp.toLowerAscii,
    &"error should mention scan, got: {outp}"
  removeDir(tmpDir)

# ---------------------------------------------------------------------------
# CL7 — testBcfScan: --scan with BCF via scatter exits 1
# ---------------------------------------------------------------------------

timed("CL7", "BCF --scan exits 1"):
  doAssert fileExists(SmallBcf), "BCF fixture missing — run generate_fixtures.sh"
  let tmpDir = createTempDir("blocky_", "")
  let (outp, code) = run(&"scatter -n 2 -o {tmpDir}/shard --scan {SmallBcf}")
  doAssert code != 0, "--scan with BCF should exit non-zero"
  doAssert "scan" in outp.toLowerAscii,
    &"error should mention scan, got: {outp}"
  removeDir(tmpDir)

# ===========================================================================
# CL8–CL13 — Integration: scatter correctness for VCF and BCF
# ===========================================================================

# ---------------------------------------------------------------------------
# CL8 — testMissingIndex: no index file → auto-scan fallback with warning; shards produced
# ---------------------------------------------------------------------------

timed("CL8", "no-index scatter: auto-scan with warning, shards produced"):
  # With no index, scatter should warn and fall back to BGZF scan automatically.
  let tmpDir = createTempDir("blocky_", "")
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
# CL9 — testScanFlag: --scan ignores existing index; shards valid, record count matches
# ---------------------------------------------------------------------------

timed("CL9", "--scan: 4 shards, records match"):
  let tmpDir = createTempDir("blocky_", "")
  let outp_template = tmpDir / "out.vcf.gz"
  let (runOutp, runCode) = run(&"scatter -n 4 -o {outp_template} --scan {SmallVcf}")
  doAssert runCode == 0, &"--scan exited non-zero:\n{runOutp}"
  proc countAndCheckScan(path: string): int =
    let (o, code) = execCmdEx("bcftools query -f '%POS\\n' " & path & " 2>/dev/null")
    doAssert code == 0, &"bcftools rejected --scan shard: {path}"
    o.strip.countLines
  var total = 0
  for i in 1..4:
    let p = tmpDir / ("shard_" & $i & ".out.vcf.gz")
    doAssert fileExists(p), &"--scan shard {i} missing"
    total += countAndCheckScan(p)
  doAssert total == countAndCheckScan(SmallVcf),
    &"--scan record count mismatch: shards={total}"
  removeDir(tmpDir)

# ---------------------------------------------------------------------------
# CL10 — testEndToEnd: scatter -n 4 VCF; all shards valid, content hash matches
# ---------------------------------------------------------------------------

timed("CL10", "e2e: record content hash matches original"):
  let tmpDir = createTempDir("blocky_", "")
  let outp_template = tmpDir / "out.vcf.gz"

  let (runOutp, runCode) = run(&"scatter -n 4 -o {outp_template} {SmallVcf}")
  doAssert runCode == 0, &"blocky scatter exited non-zero:\n{runOutp}"

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
  let tmpDir = createTempDir("blocky_", "")
  let outp_template = tmpDir / "out.vcf.gz"

  let (runOutp, runCode) = run(&"scatter -n 4 -o {outp_template} {CsiVcf}")
  doAssert runCode == 0, &"blocky scatter (CSI) exited non-zero:\n{runOutp}"

  proc countRec(path: string): int =
    let (o, _) = execCmdEx("bcftools query -f '%POS\\n' " & path & " 2>/dev/null | wc -l")
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
  let tmpDir = createTempDir("blocky_", "")
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
      let tmpDir = createTempDir("blocky_", "")
      let outp_template = tmpDir / "out.vcf.gz"

      let (runOutp, runCode) = run(&"scatter -n 10 -o {outp_template} {KgVcf}")
      doAssert runCode == 0, &"blocky scatter (1KG) exited non-zero:\n{runOutp}"

      # With -n 10, nDigits=2 so names are shard_01.out.vcf.gz … shard_10.out.vcf.gz
      for i in 1..10:
        let path = tmpDir / ("shard_" & align($i, 2, '0') & ".out.vcf.gz")
        doAssert fileExists(path), &"1KG shard {i} not found: {path}"
        let (bcfOutp, bcfCode) = execCmdEx(
          "bcftools view -HG " & path & " > /dev/null 2>&1")
        doAssert bcfCode == 0, &"bcftools rejected 1KG shard {i}: {bcfOutp}"

      proc countRecs(path: string): int =
        let (o, _) = execCmdEx("bcftools query -f '%POS\\n' " & path & " 2>/dev/null | wc -l")
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

# ---------------------------------------------------------------------------
# CL30.2 — run -u with --discard: exits non-zero (mutually exclusive)
# ---------------------------------------------------------------------------

timed("CL30.2", "run -u with --discard: exits non-zero"):
  let (outp, code) = run(
    &"run -n 2 -u --discard {SmallVcf} ::: cat")
  doAssert code != 0, &"CL30.2: -u with --discard should exit non-zero, got {code}"
  doAssert "discard" in outp.toLowerAscii,
    &"CL30.2: expected --discard error, got:\n{outp}"

# ===========================================================================
# CL31 — -n exceeds index entry count
# ===========================================================================
# Scatter and run paths must reject nShards > available index entries with a
# clear error message rather than producing empty shards. Tests cover:
#   CL31.1 — scatter VCF: error suggests --scan
#   CL31.2 — scatter BCF: error does NOT suggest --scan (BCF can't scan)
#   CL31.3 — scatter VCF --scan: succeeds when raw blocks > requested n
#   CL31.4 — run +merge+ VCF: same error path as scatter
#   CL31.5 — -n at the exact limit succeeds (boundary check)

# ---------------------------------------------------------------------------
# CL31.1 — scatter -n way too high on VCF: errors, suggests --scan
# ---------------------------------------------------------------------------

timed("CL31.1", "scatter -n too high VCF: errors with --scan suggestion"):
  let tmpDir = createTempDir("blocky_", "")
  # chr22_1kg.vcf.gz has only 57 index entries; -n 100 must error.
  let (outp, code) = run(&"scatter -n 100 -o {tmpDir}/shard_{{}}.vcf.gz {KgVcf}")
  doAssert code != 0,
    &"CL31.1: -n 100 on chr22_1kg.vcf.gz (57 voffs) should error, got {code}:\n{outp}"
  doAssert "index entries" in outp,
    &"CL31.1: error should mention 'index entries', got: {outp}"
  doAssert "--scan" in outp,
    &"CL31.1: VCF error should suggest --scan, got: {outp}"
  removeDir(tmpDir)

# ---------------------------------------------------------------------------
# CL31.2 — scatter -n way too high on BCF: errors, NO --scan suggestion
# ---------------------------------------------------------------------------

timed("CL31.2", "scatter -n too high BCF: errors without --scan suggestion"):
  doAssert fileExists(DataDir / "chr22_1kg.bcf"),
    "chr22_1kg.bcf fixture missing — run generate_fixtures.sh"
  let tmpDir = createTempDir("blocky_", "")
  # chr22_1kg.bcf has only 58 index entries; -n 100 must error.
  let (outp, code) = run(
    &"scatter -n 100 -o {tmpDir}/shard_{{}}.bcf {DataDir}/chr22_1kg.bcf")
  doAssert code != 0,
    &"CL31.2: -n 100 on chr22_1kg.bcf (58 voffs) should error, got {code}:\n{outp}"
  doAssert "index entries" in outp,
    &"CL31.2: error should mention 'index entries', got: {outp}"
  doAssert "--scan" notin outp,
    &"CL31.2: BCF error should NOT suggest --scan, got: {outp}"
  removeDir(tmpDir)

# ---------------------------------------------------------------------------
# CL31.3 — scatter --scan VCF: succeeds when raw blocks exceed -n
# ---------------------------------------------------------------------------

timed("CL31.3", "scatter --scan recovers: 30 shards, records match"):
  let tmpDir = createTempDir("blocky_", "")
  # chr22_1kg.vcf.gz has 28 unique block offsets; -n 30 exceeds index entries
  # but --scan uses raw BGZF blocks (~834 blocks) so -n 30 succeeds.
  let (outp, code) = run(
    &"scatter -n 30 --scan -o {tmpDir}/shard_{{}}.vcf.gz {KgVcf}")
  doAssert code == 0,
    &"CL31.3: --scan with -n 30 should succeed, got {code}:\n{outp}"
  let shards = toSeq(walkFiles(tmpDir / "shard_*.vcf.gz"))
  doAssert shards.len == 30,
    &"CL31.3: expected 30 shards, got {shards.len}"
  # Verify record count: concat all shards naively, count once (fast).
  let origCnt = execCmdEx(
    "bcftools query -f '%POS\\n' " & KgVcf & " 2>/dev/null | wc -l")[0].strip.parseInt
  let outCnt = execCmdEx(
    "bcftools concat --naive " & shards.join(" ") &
    " 2>/dev/null | bcftools query -f '%POS\\n' 2>/dev/null | wc -l")[0].strip.parseInt
  doAssert outCnt == origCnt,
    &"CL31.3: record count mismatch: orig={origCnt} shards={outCnt}"
  removeDir(tmpDir)

# ---------------------------------------------------------------------------
# CL31.5 — -n at the exact voff count succeeds
# ---------------------------------------------------------------------------

timed("CL31.5", "scatter -n at exact voff limit: 28 shards"):
  let tmpDir = createTempDir("blocky_", "")
  # chr22_1kg.vcf.gz has 28 voffs; -n 28 must succeed.
  let (outp, code) = run(&"scatter -n 28 -o {tmpDir}/shard_{{}}.vcf.gz {KgVcf}")
  doAssert code == 0,
    &"CL31.5: -n 28 (== voff count) should succeed, got {code}:\n{outp}"
  let shards = toSeq(walkFiles(tmpDir / "shard_*.vcf.gz"))
  doAssert shards.len == 28,
    &"CL31.5: expected 28 shards, got {shards.len}"
  removeDir(tmpDir)

# ===========================================================================
# CL32 — --clamp reduces -n instead of erroring
# ===========================================================================
# When --clamp is set and -n exceeds the available index entries, the
# tool prints an info line and silently produces fewer shards rather than
# erroring out (CL31's no-clamp behaviour).

# ---------------------------------------------------------------------------
# CL32.1 — scatter --clamp on VCF: -n 100 → 28 shards
# ---------------------------------------------------------------------------

timed("CL32.1", "scatter --clamp VCF: 28 shards, records match"):
  let tmpDir = createTempDir("blocky_", "")
  # chr22_1kg.vcf.gz has 28 index entries; -n 100 must clamp to 28.
  let (outp, code) = run(
    &"scatter -n 100 --clamp -o {tmpDir}/shard_{{}}.vcf.gz {KgVcf}")
  doAssert code == 0,
    &"CL32.1: --clamp should succeed, got {code}:\n{outp}"
  doAssert "clamp" in outp,
    &"CL32.1: info should mention 'clamp', got: {outp}"
  let shards = toSeq(walkFiles(tmpDir / "shard_*.vcf.gz"))
  doAssert shards.len == 28,
    &"CL32.1: expected 28 clamped shards, got {shards.len}"
  # Sanity-check record count matches the original.
  let origCnt = execCmdEx(
    "bcftools query -f '%POS\\n' " & KgVcf & " 2>/dev/null | wc -l")[0].strip.parseInt
  var outCnt = 0
  for s in shards:
    outCnt += execCmdEx(
      "bcftools query -f '%POS\\n' " & s & " 2>/dev/null | wc -l")[0].strip.parseInt
  doAssert outCnt == origCnt,
    &"CL32.1: record count mismatch: orig={origCnt} shards={outCnt}"
  removeDir(tmpDir)

# ===========================================================================
# CL40 — compress / decompress subcommands
# ===========================================================================

# ---------------------------------------------------------------------------
# CL40.1 — compress file, then decompress: round-trip identity
# ---------------------------------------------------------------------------

timed("CL40.1", "compress/decompress round-trip: file identity"):
  let tmpDir = createTempDir("blocky_", "")
  let rawPath = tmpDir / "small.vcf"
  # Create a raw VCF by decompressing the fixture.
  let (_, dec) = execCmdEx(&"bgzip -d -c {SmallVcf} > {rawPath}")
  doAssert dec == 0, "failed to decompress fixture"
  let origData = readFile(rawPath)

  # Compress — should create .gz and remove original.
  let (cOut, cCode) = run(&"compress {rawPath}")
  doAssert cCode == 0, &"CL40.1: compress failed: {cOut}"
  let gzPath = rawPath & ".gz"
  doAssert fileExists(gzPath), "CL40.1: compressed file not created"
  doAssert not fileExists(rawPath), "CL40.1: original not removed after compress"

  # Verify BGZF magic.
  let gzData = readFile(gzPath)
  doAssert gzData.len > 28, "CL40.1: compressed file too small"
  doAssert gzData[0].byte == 0x1f and gzData[1].byte == 0x8b,
    "CL40.1: missing gzip magic"

  # Decompress — should create raw file and remove .gz.
  let (dOut, dCode) = run(&"decompress {gzPath}")
  doAssert dCode == 0, &"CL40.1: decompress failed: {dOut}"
  doAssert fileExists(rawPath), "CL40.1: decompressed file not created"
  doAssert not fileExists(gzPath), "CL40.1: .gz not removed after decompress"
  let roundTrip = readFile(rawPath)
  doAssert roundTrip == origData,
    &"CL40.1: round-trip mismatch: orig={origData.len} got={roundTrip.len}"
  removeDir(tmpDir)

# ---------------------------------------------------------------------------
# CL40.2 — compress -c / decompress -c: stdin/stdout pipe round-trip
# ---------------------------------------------------------------------------

timed("CL40.2", "compress/decompress -c: pipe round-trip"):
  let tmpDir = createTempDir("blocky_", "")
  let rawPath = tmpDir / "small.vcf"
  discard execCmdEx(&"bgzip -d -c {SmallVcf} > {rawPath}")
  let origMd5 = execCmdEx(&"md5sum {rawPath}")[0].split(' ')[0]

  let pipeMd5 = execCmdEx(
    &"cat {rawPath} | {BinPath} compress -c | {BinPath} decompress -c | md5sum")[0].split(' ')[0]
  doAssert pipeMd5 == origMd5,
    &"CL40.2: pipe round-trip md5 mismatch: orig={origMd5} got={pipeMd5}"
  removeDir(tmpDir)

# ---------------------------------------------------------------------------
# CL40.3 — compress warns on already-compressed input
# ---------------------------------------------------------------------------

timed("CL40.3", "compress: warns on already-compressed input"):
  let (outp, code) = run(&"compress -c {SmallVcf}")
  # Should succeed but warn.
  doAssert code == 0, &"CL40.3: expected exit 0, got {code}"
  doAssert "already" in outp.toLowerAscii,
    &"CL40.3: expected warning about already-compressed, got: {outp}"

# ---------------------------------------------------------------------------
# CL40.4 — decompress warns on non-compressed input
# ---------------------------------------------------------------------------

timed("CL40.4", "decompress: warns on non-compressed input"):
  let tmpDir = createTempDir("blocky_", "")
  let rawPath = tmpDir / "small.vcf.gz"  # named .gz but actually raw
  writeFile(rawPath, "##fileformat=VCFv4.2\n#CHROM\tPOS\n")
  let (outp, code) = run(&"decompress -c {rawPath}")
  # Should warn about non-compressed input.
  doAssert "not appear" in outp.toLowerAscii or "warning" in outp.toLowerAscii,
    &"CL40.4: expected warning, got: {outp}"
  removeDir(tmpDir)

# ---------------------------------------------------------------------------
# CL40.5 — decompress: non-.gz file without -c exits 1
# ---------------------------------------------------------------------------

timed("CL40.5", "decompress: non-.gz file without -c exits 1"):
  let tmpDir = createTempDir("blocky_", "")
  let rawPath = tmpDir / "data.vcf"
  writeFile(rawPath, "test")
  let (outp, code) = run(&"decompress {rawPath}")
  doAssert code != 0, &"CL40.5: expected non-zero exit, got {code}: {outp}"
  doAssert ".gz" in outp, &"CL40.5: error should mention .gz, got: {outp}"
  removeDir(tmpDir)

# ---------------------------------------------------------------------------
# CL40.6 — compress: output already exists -> exits 1, original untouched
# ---------------------------------------------------------------------------

timed("CL40.6", "compress: output exists -> exits 1"):
  let tmpDir = createTempDir("blocky_", "")
  let rawPath = tmpDir / "small.vcf"
  let gzPath = rawPath & ".gz"
  discard execCmdEx(&"bgzip -d -c {SmallVcf} > {rawPath}")
  writeFile(gzPath, "existing")  # pre-existing output
  let origSize = getFileSize(rawPath)
  let (outp, code) = run(&"compress {rawPath}")
  doAssert code != 0, &"CL40.6: expected non-zero exit, got {code}: {outp}"
  doAssert "already exists" in outp, &"CL40.6: error should mention 'already exists', got: {outp}"
  doAssert fileExists(rawPath), "CL40.6: original should not be removed"
  doAssert getFileSize(rawPath) == origSize, "CL40.6: original should be untouched"
  doAssert readFile(gzPath) == "existing", "CL40.6: existing .gz should be untouched"
  removeDir(tmpDir)

# ---------------------------------------------------------------------------
# CL40.7 — decompress: output already exists -> exits 1, original untouched
# ---------------------------------------------------------------------------

timed("CL40.7", "decompress: output exists -> exits 1"):
  let tmpDir = createTempDir("blocky_", "")
  let rawPath = tmpDir / "small.vcf"
  let gzPath = rawPath & ".gz"
  # Create a valid .gz file to decompress.
  discard execCmdEx(&"bgzip -d -c {SmallVcf} > {rawPath}")
  discard execCmdEx(&"{BinPath} compress {rawPath}")
  doAssert fileExists(gzPath), "CL40.7: setup failed"
  # Pre-create the output file to trigger the conflict.
  writeFile(rawPath, "existing")
  let origGzSize = getFileSize(gzPath)
  let (outp, code) = run(&"decompress {gzPath}")
  doAssert code != 0, &"CL40.7: expected non-zero exit, got {code}: {outp}"
  doAssert "already exists" in outp, &"CL40.7: error should mention 'already exists', got: {outp}"
  doAssert fileExists(gzPath), "CL40.7: .gz should not be removed"
  doAssert getFileSize(gzPath) == origGzSize, "CL40.7: .gz should be untouched"
  doAssert readFile(rawPath) == "existing", "CL40.7: existing output should be untouched"
  removeDir(tmpDir)

# ===========================================================================
# CL50 — BED format: scatter, run, gather on non-VCF bgzipped file
# ===========================================================================

const SmallBed = DataDir / "small.bed.gz"     # BED, TBI indexed, no headers
const GtfFile  = DataDir / "gencode_5k.gtf.gz"  # GTF, TBI indexed, ## headers

# ---------------------------------------------------------------------------
# CL50.1 — scatter BED: 4 shards, line count matches
# ---------------------------------------------------------------------------

timed("CL50.1", "scatter BED: shards, line count matches"):
  doAssert fileExists(SmallBed), "BED fixture missing — run generate_fixtures.sh"
  let tmpDir = createTempDir("blocky_", "")
  let (outp, code) = run(
    &"scatter -n 4 --clamp -o {tmpDir}/shard_{{}}.bed.gz {SmallBed}")
  doAssert code == 0, &"CL50.1: scatter failed ({code}): {outp}"
  let shards = toSeq(walkFiles(tmpDir / "shard_*.bed.gz"))
  doAssert shards.len >= 1, &"CL50.1: expected at least 1 shard, got {shards.len}"
  let origCnt = execCmdEx(&"bgzip -d -c {SmallBed} | wc -l")[0].strip.parseInt
  var shardCnt = 0
  for s in shards:
    shardCnt += execCmdEx(&"bgzip -d -c {s} | wc -l")[0].strip.parseInt
  doAssert shardCnt == origCnt,
    &"CL50.1: line count mismatch: orig={origCnt} shards={shardCnt}"
  removeDir(tmpDir)

# ---------------------------------------------------------------------------
# CL50.2 — run BED: 4 workers, cat stage, output matches original
# ---------------------------------------------------------------------------

timed("CL50.2", "run BED: 4 workers, cat stage"):
  doAssert fileExists(SmallBed), "BED fixture missing"
  let tmpDir = createTempDir("blocky_", "")
  let outPath = tmpDir / "out.bed.gz"
  let (outp, code) = run(
    &"run -n 4 --clamp -o {outPath} {SmallBed} ::: cat")
  doAssert code == 0, &"CL50.2: run failed ({code}): {outp}"
  let origCnt = execCmdEx(&"bgzip -d -c {SmallBed} | wc -l")[0].strip.parseInt
  let outCnt = execCmdEx(&"bgzip -d -c {outPath} | wc -l")[0].strip.parseInt
  doAssert outCnt == origCnt,
    &"CL50.2: line count mismatch: orig={origCnt} out={outCnt}"
  removeDir(tmpDir)

# ---------------------------------------------------------------------------
# CL50.3 — scatter GTF: 4 shards, line count matches (skipped if no fixture)
# ---------------------------------------------------------------------------

timed("CL50.3", "scatter GTF: 4 shards, line count matches"):
  if fileExists(GtfFile):
    let tmpDir = createTempDir("blocky_", "")
    let (outp, code) = run(
      &"scatter -n 4 --clamp -o {tmpDir}/shard_{{}}.gtf.gz {GtfFile}")
    doAssert code == 0, &"CL50.3: scatter failed ({code}): {outp}"
    let shards = toSeq(walkFiles(tmpDir / "shard_*.gtf.gz"))
    doAssert shards.len >= 1, &"CL50.3: expected at least 1 shard, got {shards.len}"
    # Count data lines only (exclude # headers — each shard includes headers).
    let origCnt = execCmdEx(&"bgzip -d -c {GtfFile} | grep -cv '^#'")[0].strip.parseInt
    var shardCnt = 0
    for s in shards:
      shardCnt += execCmdEx(&"bgzip -d -c {s} | grep -cv '^#'")[0].strip.parseInt
    doAssert shardCnt == origCnt,
      &"CL50.3: data line count mismatch: orig={origCnt} shards={shardCnt}"
    removeDir(tmpDir)
  else:
    echo "  [skipped — gencode_5k.gtf.gz not found]"

# ===========================================================================
# CL51 — --scan on unindexed file
# ===========================================================================

const SmallTsv = DataDir / "small_noindex.tsv.gz"  # unindexed, # header

# ---------------------------------------------------------------------------
# CL51.1 — scatter --scan unindexed TSV: 4 shards, line count matches
# ---------------------------------------------------------------------------

timed("CL51.1", "scatter --scan unindexed TSV: line count matches"):
  doAssert fileExists(SmallTsv), "TSV fixture missing — run generate_fixtures.sh"
  let tmpDir = createTempDir("blocky_", "")
  let (outp, code) = run(
    &"scatter -n 4 --scan --clamp -o {tmpDir}/shard_{{}}.tsv.gz {SmallTsv}")
  doAssert code == 0, &"CL51.1: scatter failed ({code}): {outp}"
  let shards = toSeq(walkFiles(tmpDir / "shard_*.tsv.gz"))
  doAssert shards.len >= 1, &"CL51.1: expected at least 1 shard, got {shards.len}"
  # Count data lines only (exclude # headers — each shard includes headers).
  let origCnt = execCmdEx(&"bgzip -d -c {SmallTsv} | grep -cv '^#'")[0].strip.parseInt
  var shardCnt = 0
  for s in shards:
    shardCnt += execCmdEx(&"bgzip -d -c {s} | grep -cv '^#'")[0].strip.parseInt
  doAssert shardCnt == origCnt,
    &"CL51.1: data line count mismatch: orig={origCnt} shards={shardCnt}"
  removeDir(tmpDir)

# ---------------------------------------------------------------------------
# CL51.2 — auto-scan unindexed TSV: warning + shards produced
# ---------------------------------------------------------------------------

timed("CL51.2", "auto-scan unindexed TSV: warning + shards produced"):
  doAssert fileExists(SmallTsv), "TSV fixture missing"
  let tmpDir = createTempDir("blocky_", "")
  let (outp, code) = run(
    &"scatter -n 4 --clamp -o {tmpDir}/shard_{{}}.tsv.gz {SmallTsv}")
  doAssert code == 0, &"CL51.2: scatter failed ({code}): {outp}"
  doAssert "warning" in outp.toLowerAscii,
    &"CL51.2: expected warning about missing index, got: {outp}"
  let shards = toSeq(walkFiles(tmpDir / "shard_*.tsv.gz"))
  doAssert shards.len >= 1, &"CL51.2: expected at least 1 shard, got {shards.len}"
  removeDir(tmpDir)

