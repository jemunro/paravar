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

timed("C01", "binary available"):
  if not fileExists(BinPath):
    let (outp, code) = execCmdEx("nimble build 2>&1")
    if code != 0:
      echo "nimble build failed:\n", outp
      quit(1)
  doAssert fileExists(BinPath), "binary not found: " & BinPath & " (run nimble build)"

# ===========================================================================
# C02–C08 — Error cases (missing flags, invalid args, unsupported modes)
# ===========================================================================

# ---------------------------------------------------------------------------
# C02 — testMissingN: missing -n exits non-zero with 'n' in error message
# ---------------------------------------------------------------------------

timed("C02", "missing -n exits non-zero"):
  let (outp, code) = run(&"scatter -o /tmp/blocky_cli_test {SmallVcf}")
  doAssert code != 0, "missing -n should exit non-zero"
  doAssert "n" in outp.toLowerAscii,
    &"missing -n error should mention 'n', got: {outp}"

# ---------------------------------------------------------------------------
# C03 — testInvalidN0: -n 0 exits non-zero
# ---------------------------------------------------------------------------

timed("C03", "-n 0 exits non-zero"):
  let (_, code) = run(&"scatter -n 0 -o /tmp/blocky_cli_test {SmallVcf}")
  doAssert code != 0, "-n 0 should exit non-zero"

# ---------------------------------------------------------------------------
# C04 — testMissingO: missing -o exits non-zero with 'o' in error message
# ---------------------------------------------------------------------------

timed("C04", "missing -o exits non-zero"):
  let (outp, code) = run(&"scatter -n 2 {SmallVcf}")
  doAssert code != 0, "missing -o should exit non-zero"
  doAssert "o" in outp.toLowerAscii,
    &"missing -o error should mention 'o', got: {outp}"

# ---------------------------------------------------------------------------
# C05 — testUnknownExtension: unknown extension (.xyz) exits 1 with extension in message
# ---------------------------------------------------------------------------

timed("C05", "unknown extension exits 1 with extension in message"):
  let tmpDir = createTempDir("blocky_", "")
  let tmpFile = tmpDir / "input.xyz"
  writeFile(tmpFile, "dummy")
  let (outp, code) = run(&"scatter -n 2 -o {tmpDir}/shard {tmpFile}")
  doAssert code != 0, "unknown extension should exit non-zero"
  doAssert ".xyz" in outp, &"error message should contain '.xyz', got: {outp}"
  removeDir(tmpDir)

# ---------------------------------------------------------------------------
# C06 — testBcfNoIndex: BCF without index exits 1 (no auto-scan fallback)
# ---------------------------------------------------------------------------

timed("C06", "BCF no index exits 1"):
  doAssert fileExists(SmallBcf), "BCF fixture missing — run generate_fixtures.sh"
  let tmpDir = createTempDir("blocky_", "")
  let tmpBcf = tmpDir / "noindex.bcf"
  copyFile(SmallBcf, tmpBcf)
  # Intentionally omit .csi so BCF has no index.
  let (outp, code) = run(&"scatter -n 2 -o {tmpDir}/shard {tmpBcf}")
  doAssert code != 0, &"BCF with no index should exit non-zero, got {code}:\n{outp}"
  removeDir(tmpDir)

# ---------------------------------------------------------------------------
# C07 — testBcfRunScan: --scan with BCF via run exits 1
# ---------------------------------------------------------------------------

timed("C07", "BCF run --scan exits 1"):
  doAssert fileExists(SmallBcf), "BCF fixture missing — run generate_fixtures.sh"
  let tmpDir = createTempDir("blocky_", "")
  let (outp, code) = run(&"run -n 2 -o {tmpDir}/out.vcf.gz --scan {SmallBcf} --- cat")
  doAssert code != 0, "--scan with BCF via run should exit non-zero"
  doAssert "scan" in outp.toLowerAscii,
    &"error should mention scan, got: {outp}"
  removeDir(tmpDir)

# ---------------------------------------------------------------------------
# C08 — testBcfScan: --scan with BCF via scatter exits 1
# ---------------------------------------------------------------------------

timed("C08", "BCF --scan exits 1"):
  doAssert fileExists(SmallBcf), "BCF fixture missing — run generate_fixtures.sh"
  let tmpDir = createTempDir("blocky_", "")
  let (outp, code) = run(&"scatter -n 2 -o {tmpDir}/shard --scan {SmallBcf}")
  doAssert code != 0, "--scan with BCF should exit non-zero"
  doAssert "scan" in outp.toLowerAscii,
    &"error should mention scan, got: {outp}"
  removeDir(tmpDir)

# ===========================================================================
# C09–C14 — Integration: scatter correctness for VCF and BCF
# ===========================================================================

# ---------------------------------------------------------------------------
# C09 — testMissingIndex: no index file → auto-scan fallback with warning; shards produced
# ---------------------------------------------------------------------------

timed("C09", "no-index scatter: auto-scan with warning, shards produced"):
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
# C10 — testScanFlag: --scan ignores existing index; shards valid, record count matches
# ---------------------------------------------------------------------------

timed("C10", "--scan: 4 shards, records match"):
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
# C11 — testEndToEnd: scatter -n 4 VCF; all shards valid, content hash matches
# ---------------------------------------------------------------------------

timed("C11", "e2e: record content hash matches original"):
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
# C12 — testCsiIndex: CSI-indexed VCF scatter; shards valid, record count matches
# ---------------------------------------------------------------------------

timed("C12", "CSI: scatter -n 4, all shards valid"):
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
# C13 — testBcfExtension: BCF scatter produces .bcf shards with matching content hash
# ---------------------------------------------------------------------------

timed("C13", "BCF .bcf extension -> BCF mode, content hash matches"):
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
# C14 — testKg1000Genomes: large 1KG VCF scatter -n 10 (skipped if fixture absent)
# ---------------------------------------------------------------------------

block testKg1000Genomes:
  if not fileExists(KgVcf):
    echo "SKIP C14 1KG chr22: file not present (run tests/generate_fixtures.sh)"
  else:
    timed("C14", "1KG chr22 scatter: records match original"):
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
# C15 — -u flag (force uncompressed output)
# ===========================================================================

# ---------------------------------------------------------------------------
# C15 — run -u with --discard: exits non-zero (mutually exclusive)
# ---------------------------------------------------------------------------

timed("C15", "run -u with --discard: exits non-zero"):
  let (outp, code) = run(
    &"run -n 2 -u --discard {SmallVcf} ::: cat")
  doAssert code != 0, &"C15: -u with --discard should exit non-zero, got {code}"
  doAssert "discard" in outp.toLowerAscii,
    &"C15: expected --discard error, got:\n{outp}"

# ===========================================================================
# C025 — -n exceeds index entry count
# ===========================================================================
# Scatter and run paths must reject nShards > available index entries with a
# clear error message rather than producing empty shards. Tests cover:
#   C16 — scatter VCF: error suggests --scan
#   C17 — scatter BCF: error does NOT suggest --scan (BCF can't scan)
#   C18 — scatter VCF --scan: succeeds when raw blocks > requested n
#   C19 — -n at the exact limit succeeds (boundary check)

# ---------------------------------------------------------------------------
# C16 — scatter -n way too high on VCF: errors, suggests --scan
# ---------------------------------------------------------------------------

timed("C16", "scatter -n too high VCF: errors with --scan suggestion"):
  let tmpDir = createTempDir("blocky_", "")
  # chr22_1kg.vcf.gz has only 57 index entries; -n 100 must error.
  let (outp, code) = run(&"scatter -n 100 -o {tmpDir}/shard_{{}}.vcf.gz {KgVcf}")
  doAssert code != 0,
    &"C16: -n 100 on chr22_1kg.vcf.gz (57 voffs) should error, got {code}:\n{outp}"
  doAssert "index entries" in outp,
    &"C16: error should mention 'index entries', got: {outp}"
  doAssert "--scan" in outp,
    &"C16: VCF error should suggest --scan, got: {outp}"
  removeDir(tmpDir)

# ---------------------------------------------------------------------------
# C17 — scatter -n way too high on BCF: errors, NO --scan suggestion
# ---------------------------------------------------------------------------

timed("C17", "scatter -n too high BCF: errors without --scan suggestion"):
  doAssert fileExists(DataDir / "chr22_1kg.bcf"),
    "chr22_1kg.bcf fixture missing — run generate_fixtures.sh"
  let tmpDir = createTempDir("blocky_", "")
  # chr22_1kg.bcf has only 58 index entries; -n 100 must error.
  let (outp, code) = run(
    &"scatter -n 100 -o {tmpDir}/shard_{{}}.bcf {DataDir}/chr22_1kg.bcf")
  doAssert code != 0,
    &"C17: -n 100 on chr22_1kg.bcf (58 voffs) should error, got {code}:\n{outp}"
  doAssert "index entries" in outp,
    &"C17: error should mention 'index entries', got: {outp}"
  doAssert "--scan" notin outp,
    &"C17: BCF error should NOT suggest --scan, got: {outp}"
  removeDir(tmpDir)

# ---------------------------------------------------------------------------
# C18 — scatter --scan VCF: succeeds when raw blocks exceed -n
# ---------------------------------------------------------------------------

timed("C18", "scatter --scan recovers: 30 shards, records match"):
  let tmpDir = createTempDir("blocky_", "")
  # chr22_1kg.vcf.gz has 28 unique block offsets; -n 30 exceeds index entries
  # but --scan uses raw BGZF blocks (~834 blocks) so -n 30 succeeds.
  let (outp, code) = run(
    &"scatter -n 30 --scan -o {tmpDir}/shard_{{}}.vcf.gz {KgVcf}")
  doAssert code == 0,
    &"C18: --scan with -n 30 should succeed, got {code}:\n{outp}"
  let shards = toSeq(walkFiles(tmpDir / "shard_*.vcf.gz"))
  doAssert shards.len == 30,
    &"C18: expected 30 shards, got {shards.len}"
  # Verify record count: concat all shards naively, count once (fast).
  let origCnt = execCmdEx(
    "bcftools query -f '%POS\\n' " & KgVcf & " 2>/dev/null | wc -l")[0].strip.parseInt
  let outCnt = execCmdEx(
    "bcftools concat --naive " & shards.join(" ") &
    " 2>/dev/null | bcftools query -f '%POS\\n' 2>/dev/null | wc -l")[0].strip.parseInt
  doAssert outCnt == origCnt,
    &"C18: record count mismatch: orig={origCnt} shards={outCnt}"
  removeDir(tmpDir)

# ---------------------------------------------------------------------------
# C19 — -n at the exact voff count succeeds
# ---------------------------------------------------------------------------

timed("C19", "scatter -n at exact voff limit: 28 shards"):
  let tmpDir = createTempDir("blocky_", "")
  # chr22_1kg.vcf.gz has 28 voffs; -n 28 must succeed.
  let (outp, code) = run(&"scatter -n 28 -o {tmpDir}/shard_{{}}.vcf.gz {KgVcf}")
  doAssert code == 0,
    &"C19: -n 28 (== voff count) should succeed, got {code}:\n{outp}"
  let shards = toSeq(walkFiles(tmpDir / "shard_*.vcf.gz"))
  doAssert shards.len == 28,
    &"C19: expected 28 shards, got {shards.len}"
  removeDir(tmpDir)

# ===========================================================================
# C20 — --clamp reduces -n instead of erroring
# ===========================================================================
# When --clamp is set and -n exceeds the available index entries, the
# tool prints an info line and silently produces fewer shards rather than
# erroring out (C025's no-clamp behaviour).

# ---------------------------------------------------------------------------
# C20 — scatter --clamp on VCF: -n 100 → 28 shards
# ---------------------------------------------------------------------------

timed("C20", "scatter --clamp VCF: 28 shards, records match"):
  let tmpDir = createTempDir("blocky_", "")
  # chr22_1kg.vcf.gz has 28 index entries; -n 100 must clamp to 28.
  let (outp, code) = run(
    &"scatter -n 100 --clamp -o {tmpDir}/shard_{{}}.vcf.gz {KgVcf}")
  doAssert code == 0,
    &"C20: --clamp should succeed, got {code}:\n{outp}"
  doAssert "clamp" in outp,
    &"C20: info should mention 'clamp', got: {outp}"
  let shards = toSeq(walkFiles(tmpDir / "shard_*.vcf.gz"))
  doAssert shards.len == 28,
    &"C20: expected 28 clamped shards, got {shards.len}"
  # Sanity-check record count matches the original.
  let origCnt = execCmdEx(
    "bcftools query -f '%POS\\n' " & KgVcf & " 2>/dev/null | wc -l")[0].strip.parseInt
  var outCnt = 0
  for s in shards:
    outCnt += execCmdEx(
      "bcftools query -f '%POS\\n' " & s & " 2>/dev/null | wc -l")[0].strip.parseInt
  doAssert outCnt == origCnt,
    &"C20: record count mismatch: orig={origCnt} shards={outCnt}"
  removeDir(tmpDir)

# ===========================================================================
# C027 — compress / decompress subcommands
# ===========================================================================

# ---------------------------------------------------------------------------
# C21 — compress file, then decompress: round-trip identity
# ---------------------------------------------------------------------------

timed("C21", "compress/decompress round-trip: file identity"):
  let tmpDir = createTempDir("blocky_", "")
  let rawPath = tmpDir / "small.vcf"
  # Create a raw VCF by decompressing the fixture.
  let (_, dec) = execCmdEx(&"bgzip -d -c {SmallVcf} > {rawPath}")
  doAssert dec == 0, "failed to decompress fixture"
  let origData = readFile(rawPath)

  # Compress — should create .gz and remove original.
  let (cOut, cCode) = run(&"compress {rawPath}")
  doAssert cCode == 0, &"C21: compress failed: {cOut}"
  let gzPath = rawPath & ".gz"
  doAssert fileExists(gzPath), "C21: compressed file not created"
  doAssert not fileExists(rawPath), "C21: original not removed after compress"

  # Verify BGZF magic.
  let gzData = readFile(gzPath)
  doAssert gzData.len > 28, "C21: compressed file too small"
  doAssert gzData[0].byte == 0x1f and gzData[1].byte == 0x8b,
    "C21: missing gzip magic"

  # Decompress — should create raw file and remove .gz.
  let (dOut, dCode) = run(&"decompress {gzPath}")
  doAssert dCode == 0, &"C21: decompress failed: {dOut}"
  doAssert fileExists(rawPath), "C21: decompressed file not created"
  doAssert not fileExists(gzPath), "C21: .gz not removed after decompress"
  let roundTrip = readFile(rawPath)
  doAssert roundTrip == origData,
    &"C21: round-trip mismatch: orig={origData.len} got={roundTrip.len}"
  removeDir(tmpDir)

# ---------------------------------------------------------------------------
# C22 — compress -c / decompress -c: stdin/stdout pipe round-trip
# ---------------------------------------------------------------------------

timed("C22", "compress/decompress -c: pipe round-trip"):
  let tmpDir = createTempDir("blocky_", "")
  let rawPath = tmpDir / "small.vcf"
  discard execCmdEx(&"bgzip -d -c {SmallVcf} > {rawPath}")
  let origMd5 = execCmdEx(&"md5sum {rawPath}")[0].split(' ')[0]

  let pipeMd5 = execCmdEx(
    &"cat {rawPath} | {BinPath} compress -c | {BinPath} decompress -c | md5sum")[0].split(' ')[0]
  doAssert pipeMd5 == origMd5,
    &"C22: pipe round-trip md5 mismatch: orig={origMd5} got={pipeMd5}"
  removeDir(tmpDir)

# ---------------------------------------------------------------------------
# C23 — compress warns on already-compressed input
# ---------------------------------------------------------------------------

timed("C23", "compress: warns on already-compressed input"):
  let (outp, code) = run(&"compress -c {SmallVcf}")
  # Should succeed but warn.
  doAssert code == 0, &"C23: expected exit 0, got {code}"
  doAssert "already" in outp.toLowerAscii,
    &"C23: expected warning about already-compressed, got: {outp}"

# ---------------------------------------------------------------------------
# C24 — decompress warns on non-compressed input
# ---------------------------------------------------------------------------

timed("C24", "decompress: warns on non-compressed input"):
  let tmpDir = createTempDir("blocky_", "")
  let rawPath = tmpDir / "small.vcf.gz"  # named .gz but actually raw
  writeFile(rawPath, "##fileformat=VCFv4.2\n#CHROM\tPOS\n")
  let (outp, code) = run(&"decompress -c {rawPath}")
  # Should warn about non-compressed input.
  doAssert "not appear" in outp.toLowerAscii or "warning" in outp.toLowerAscii,
    &"C24: expected warning, got: {outp}"
  removeDir(tmpDir)

# ---------------------------------------------------------------------------
# C25 — decompress: non-.gz file without -c exits 1
# ---------------------------------------------------------------------------

timed("C25", "decompress: non-.gz file without -c exits 1"):
  let tmpDir = createTempDir("blocky_", "")
  let rawPath = tmpDir / "data.vcf"
  writeFile(rawPath, "test")
  let (outp, code) = run(&"decompress {rawPath}")
  doAssert code != 0, &"C25: expected non-zero exit, got {code}: {outp}"
  doAssert ".gz" in outp, &"C25: error should mention .gz, got: {outp}"
  removeDir(tmpDir)

# ---------------------------------------------------------------------------
# C26 — compress: output already exists -> exits 1, original untouched
# ---------------------------------------------------------------------------

timed("C26", "compress: output exists -> exits 1"):
  let tmpDir = createTempDir("blocky_", "")
  let rawPath = tmpDir / "small.vcf"
  let gzPath = rawPath & ".gz"
  discard execCmdEx(&"bgzip -d -c {SmallVcf} > {rawPath}")
  writeFile(gzPath, "existing")  # pre-existing output
  let origSize = getFileSize(rawPath)
  let (outp, code) = run(&"compress {rawPath}")
  doAssert code != 0, &"C26: expected non-zero exit, got {code}: {outp}"
  doAssert "already exists" in outp, &"C26: error should mention 'already exists', got: {outp}"
  doAssert fileExists(rawPath), "C26: original should not be removed"
  doAssert getFileSize(rawPath) == origSize, "C26: original should be untouched"
  doAssert readFile(gzPath) == "existing", "C26: existing .gz should be untouched"
  removeDir(tmpDir)

# ---------------------------------------------------------------------------
# C27 — decompress: output already exists -> exits 1, original untouched
# ---------------------------------------------------------------------------

timed("C27", "decompress: output exists -> exits 1"):
  let tmpDir = createTempDir("blocky_", "")
  let rawPath = tmpDir / "small.vcf"
  let gzPath = rawPath & ".gz"
  # Create a valid .gz file to decompress.
  discard execCmdEx(&"bgzip -d -c {SmallVcf} > {rawPath}")
  discard execCmdEx(&"{BinPath} compress {rawPath}")
  doAssert fileExists(gzPath), "C27: setup failed"
  # Pre-create the output file to trigger the conflict.
  writeFile(rawPath, "existing")
  let origGzSize = getFileSize(gzPath)
  let (outp, code) = run(&"decompress {gzPath}")
  doAssert code != 0, &"C27: expected non-zero exit, got {code}: {outp}"
  doAssert "already exists" in outp, &"C27: error should mention 'already exists', got: {outp}"
  doAssert fileExists(gzPath), "C27: .gz should not be removed"
  doAssert getFileSize(gzPath) == origGzSize, "C27: .gz should be untouched"
  doAssert readFile(rawPath) == "existing", "C27: existing output should be untouched"
  removeDir(tmpDir)

# ===========================================================================
# C028 — BED format: scatter, run, gather on non-VCF bgzipped file
# ===========================================================================

const SmallBed = DataDir / "small.bed.gz"     # BED, TBI indexed, no headers
const GtfFile  = DataDir / "gencode_5k.gtf.gz"  # GTF, TBI indexed, ## headers

# ---------------------------------------------------------------------------
# C28 — scatter BED: 4 shards, line count matches
# ---------------------------------------------------------------------------

timed("C28", "scatter BED: shards, line count matches"):
  doAssert fileExists(SmallBed), "BED fixture missing — run generate_fixtures.sh"
  let tmpDir = createTempDir("blocky_", "")
  let (outp, code) = run(
    &"scatter -n 4 --clamp -o {tmpDir}/shard_{{}}.bed.gz {SmallBed}")
  doAssert code == 0, &"C28: scatter failed ({code}): {outp}"
  let shards = toSeq(walkFiles(tmpDir / "shard_*.bed.gz"))
  doAssert shards.len >= 1, &"C28: expected at least 1 shard, got {shards.len}"
  let origCnt = execCmdEx(&"bgzip -d -c {SmallBed} | wc -l")[0].strip.parseInt
  var shardCnt = 0
  for s in shards:
    shardCnt += execCmdEx(&"bgzip -d -c {s} | wc -l")[0].strip.parseInt
  doAssert shardCnt == origCnt,
    &"C28: line count mismatch: orig={origCnt} shards={shardCnt}"
  removeDir(tmpDir)

# ---------------------------------------------------------------------------
# C29 — run BED: 4 workers, cat stage, output matches original
# ---------------------------------------------------------------------------

timed("C29", "run BED: 4 workers, cat stage"):
  doAssert fileExists(SmallBed), "BED fixture missing"
  let tmpDir = createTempDir("blocky_", "")
  let outPath = tmpDir / "out.bed.gz"
  let (outp, code) = run(
    &"run -n 4 --clamp -o {outPath} {SmallBed} ::: cat")
  doAssert code == 0, &"C29: run failed ({code}): {outp}"
  let origCnt = execCmdEx(&"bgzip -d -c {SmallBed} | wc -l")[0].strip.parseInt
  let outCnt = execCmdEx(&"bgzip -d -c {outPath} | wc -l")[0].strip.parseInt
  doAssert outCnt == origCnt,
    &"C29: line count mismatch: orig={origCnt} out={outCnt}"
  removeDir(tmpDir)

# ---------------------------------------------------------------------------
# C30 — scatter GTF: 4 shards, line count matches (skipped if no fixture)
# ---------------------------------------------------------------------------

timed("C30", "scatter GTF: 4 shards, line count matches"):
  if fileExists(GtfFile):
    let tmpDir = createTempDir("blocky_", "")
    let (outp, code) = run(
      &"scatter -n 4 --clamp -o {tmpDir}/shard_{{}}.gtf.gz {GtfFile}")
    doAssert code == 0, &"C30: scatter failed ({code}): {outp}"
    let shards = toSeq(walkFiles(tmpDir / "shard_*.gtf.gz"))
    doAssert shards.len >= 1, &"C30: expected at least 1 shard, got {shards.len}"
    # Count data lines only (exclude # headers — each shard includes headers).
    let origCnt = execCmdEx(&"bgzip -d -c {GtfFile} | grep -cv '^#'")[0].strip.parseInt
    var shardCnt = 0
    for s in shards:
      shardCnt += execCmdEx(&"bgzip -d -c {s} | grep -cv '^#'")[0].strip.parseInt
    doAssert shardCnt == origCnt,
      &"C30: data line count mismatch: orig={origCnt} shards={shardCnt}"
    removeDir(tmpDir)
  else:
    echo "  [skipped — gencode_5k.gtf.gz not found]"

# ===========================================================================
# C029 — --scan on unindexed file
# ===========================================================================

const SmallTsv = DataDir / "small_noindex.tsv.gz"  # unindexed, # header

# ---------------------------------------------------------------------------
# C31 — scatter --scan unindexed TSV: 4 shards, line count matches
# ---------------------------------------------------------------------------

timed("C31", "scatter --scan unindexed TSV: line count matches"):
  doAssert fileExists(SmallTsv), "TSV fixture missing — run generate_fixtures.sh"
  let tmpDir = createTempDir("blocky_", "")
  let (outp, code) = run(
    &"scatter -n 4 --scan --clamp -o {tmpDir}/shard_{{}}.tsv.gz {SmallTsv}")
  doAssert code == 0, &"C31: scatter failed ({code}): {outp}"
  let shards = toSeq(walkFiles(tmpDir / "shard_*.tsv.gz"))
  doAssert shards.len >= 1, &"C31: expected at least 1 shard, got {shards.len}"
  # Count data lines only (exclude # headers — each shard includes headers).
  let origCnt = execCmdEx(&"bgzip -d -c {SmallTsv} | grep -cv '^#'")[0].strip.parseInt
  var shardCnt = 0
  for s in shards:
    shardCnt += execCmdEx(&"bgzip -d -c {s} | grep -cv '^#'")[0].strip.parseInt
  doAssert shardCnt == origCnt,
    &"C31: data line count mismatch: orig={origCnt} shards={shardCnt}"
  removeDir(tmpDir)

# ---------------------------------------------------------------------------
# C32 — auto-scan unindexed TSV: warning + shards produced
# ---------------------------------------------------------------------------

timed("C32", "auto-scan unindexed TSV: warning + shards produced"):
  doAssert fileExists(SmallTsv), "TSV fixture missing"
  let tmpDir = createTempDir("blocky_", "")
  let (outp, code) = run(
    &"scatter -n 4 --clamp -o {tmpDir}/shard_{{}}.tsv.gz {SmallTsv}")
  doAssert code == 0, &"C32: scatter failed ({code}): {outp}"
  doAssert "warning" in outp.toLowerAscii,
    &"C32: expected warning about missing index, got: {outp}"
  let shards = toSeq(walkFiles(tmpDir / "shard_*.tsv.gz"))
  doAssert shards.len >= 1, &"C32: expected at least 1 shard, got {shards.len}"
  removeDir(tmpDir)

