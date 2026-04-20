## Tests for gather.nim — gather subcommand, #CHROM validation.
## Run from project root: nim c -d:debug -r tests/test_gather.nim

echo "--------------- Test Gather ---------------"

import std/[os, osproc, strformat, strutils, tempfiles]
import std/posix
import test_utils
import "../src/blocky/gather"
import "../src/blocky/bgzf"

const DataDir  = "tests/data"
const SmallVcf = DataDir / "small.vcf.gz"
const SmallBcf = DataDir / "small.bcf"

# ---------------------------------------------------------------------------
# G01 — testInferAllExtensions: all valid extensions (happy paths)
# ---------------------------------------------------------------------------

timed("G01", "inferFileFormat: all extensions"):
  block:
    let (f, c) = inferFileFormat("out.vcf.gz", "")
    doAssert f == ffVcf,          &".vcf.gz: expected ffVcf, got {f}"
    doAssert c == compBgzf,         &".vcf.gz: expected compBgzf, got {c}"
  block:
    let (f, c) = inferFileFormat("out.vcf.bgz", "")
    doAssert f == ffVcf,          &".vcf.bgz: expected ffVcf, got {f}"
    doAssert c == compBgzf,         &".vcf.bgz: expected compBgzf, got {c}"
  block:
    let (f, c) = inferFileFormat("out.bcf", "")
    doAssert f == ffBcf,          &".bcf: expected ffBcf, got {f}"
    doAssert c == compBgzf,         &".bcf: expected compBgzf, got {c}"
  block:
    let (f, c) = inferFileFormat("out.vcf", "")
    doAssert f == ffVcf,          &".vcf: expected ffVcf, got {f}"
    doAssert c == compNone, &".vcf: expected compNone, got {c}"
  block:
    let (f, c) = inferFileFormat("out.txt.gz", "")
    doAssert f == ffText,         &".txt.gz: expected ffText, got {f}"
    doAssert c == compBgzf,         &".txt.gz: expected compBgzf, got {c}"
  block:
    let (f, c) = inferFileFormat("out.txt", "")
    doAssert f == ffText,         &".txt: expected ffText, got {f}"
    doAssert c == compNone, &".txt: expected compNone, got {c}"
  block:
    let (f, c) = inferFileFormat("out.out.gz", "")
    doAssert f == ffText,         &"any.gz: expected ffText, got {f}"
    doAssert c == compBgzf,         &"any.gz: expected compBgzf, got {c}"
  block:
    let (f, c) = inferFileFormat("out.bgz", "")
    doAssert f == ffText,         &".bgz: expected ffText, got {f}"
    doAssert c == compBgzf,         &".bgz: expected compBgzf, got {c}"
  block:
    let (f, c) = inferFileFormat("out.xyz", "")
    doAssert f == ffText,         &".xyz: expected ffText (unknown ext), got {f}"
    doAssert c == compNone, &".xyz: expected compNone, got {c}"

# ---------------------------------------------------------------------------
# G02 — testInferWithOverride: fmtOverride overrides format, compression from extension
# ---------------------------------------------------------------------------

timed("G02", "inferFileFormat: fmtOverride overrides format"):
  # Override format; compression still from path extension.
  block:
    let (f, c) = inferFileFormat("out.vcf.gz", "bcf")
    doAssert f == ffBcf and c == compBgzf,
      &".vcf.gz + --gather-fmt bcf: got ({f}, {c})"
  block:
    let (f, c) = inferFileFormat("out.vcf.gz", "txt")
    doAssert f == ffText and c == compBgzf,
      &".vcf.gz + --gather-fmt txt: got ({f}, {c})"
  # Unknown extension with valid override → format from override, compression compNone.
  block:
    let (f, c) = inferFileFormat("out.xyz", "vcf")
    doAssert f == ffVcf and c == compNone,
      &".xyz + --gather-fmt vcf: got ({f}, {c})"

# ---------------------------------------------------------------------------
# G03 — testFindBcfHeaderEnd: offset past magic+l_text+header text returned correctly
# ---------------------------------------------------------------------------

timed("G03", "findBcfHeaderEnd"):
  # Buffer shorter than 9 bytes → -1.
  doAssert findBcfHeaderEnd(@[0x42'u8, 0x43, 0x46]) == -1, "short: expected -1"

  # Full header: magic(5) + l_text(4) + text.
  let hdrText = "##BCF header\n"
  let lText = hdrText.len.uint32
  var buf: seq[byte]
  buf.add([byte('B'), byte('C'), byte('F'), 0x02'u8, 0x02'u8])
  buf.add([byte(lText and 0xff), byte((lText shr 8) and 0xff),
           byte((lText shr 16) and 0xff), byte((lText shr 24) and 0xff)])
  for c in hdrText: buf.add(byte(c))
  let expected = 5 + 4 + hdrText.len
  doAssert findBcfHeaderEnd(buf) == expected,
    &"full header: expected {expected}, got {findBcfHeaderEnd(buf)}"

  # Partial header (one byte short) → -1.
  doAssert findBcfHeaderEnd(buf[0 ..< buf.len - 1]) == -1,
    "partial header: expected -1"

  # Extra record bytes after header → still returns headerEnd.
  var buf2 = buf & @[0x01'u8, 0x02, 0x03]
  doAssert findBcfHeaderEnd(buf2) == expected,
    "header + records: expected same headerEnd"

# ---------------------------------------------------------------------------
# G04 — testFindVcfHeaderEnd: data start offset returned; all-hash returns -1
# ---------------------------------------------------------------------------

timed("G04", "findVcfHeaderEnd"):
  # All '#' lines, no data → -1.
  var allHash: seq[byte]
  for c in "##format=VCFv4.2\n#CHROM\tPOS\n": allHash.add(byte(c))
  doAssert findVcfHeaderEnd(allHash) == -1, "all-hash lines: expected -1"

  # '#' lines followed by a data line → returns start of data line.
  var withData: seq[byte]
  for c in "##format=VCFv4.2\n#CHROM\tPOS\n1\t100\n": withData.add(byte(c))
  let dataStart = "##format=VCFv4.2\n#CHROM\tPOS\n".len
  doAssert findVcfHeaderEnd(withData) == dataStart,
    &"with data: expected {dataStart}, got {findVcfHeaderEnd(withData)}"

  # First byte is data (no header) → returns 0.
  var noHeader: seq[byte]
  for c in "1\t100\n2\t200\n": noHeader.add(byte(c))
  doAssert findVcfHeaderEnd(noHeader) == 0, "no header: expected 0"

  # Ends mid '#' line (no newline yet) → -1.
  var midLine: seq[byte]
  for c in "##format": midLine.add(byte(c))
  doAssert findVcfHeaderEnd(midLine) == -1, "mid-hash-line: expected -1"

  # Ends mid data line (no newline yet) → returns start of that line.
  var midData: seq[byte]
  for c in "##hdr\n1\t10": midData.add(byte(c))
  let midDataStart = "##hdr\n".len
  doAssert findVcfHeaderEnd(midData) == midDataStart,
    &"mid-data-line: expected {midDataStart}, got {findVcfHeaderEnd(midData)}"

# ===========================================================================
# G05 — Binary setup
# ===========================================================================

const BinPath = "./blocky"

timed("G05", "binary available"):
  if not fileExists(BinPath):
    let (outp, code) = execCmdEx("nimble build 2>&1")
    if code != 0:
      echo "nimble build failed:\n", outp
      quit(1)
  doAssert fileExists(BinPath), "binary not found: " & BinPath

proc countRecords(path: string): int =
  let (o, _) = execCmdEx("bcftools query -f '%POS\\n' " & path & " 2>/dev/null | wc -l")
  o.strip.parseInt

proc recordsHash(path: string): string =
  let (h, _) = execCmdEx("bcftools view -H " & path & " 2>/dev/null | sha256sum")
  h.split(" ")[0]

# ===========================================================================
# G06–G10 — gather subcommand integration tests
# ===========================================================================

proc runGatherSubcmd(args: string): (string, int) =
  execCmdEx(BinPath & " gather " & args & " 2>&1")

# ---------------------------------------------------------------------------
# G06 — VCF: scatter → gather, record count and hash match
# ---------------------------------------------------------------------------

timed("G06", "gather subcommand VCF: records, content hash matches"):
  let tmpDir = createTempDir("blocky_", "")
  let shardTemplate = tmpDir / "shard.{}.vcf.gz"
  # Scatter into 4 shards via CLI.
  let (sOutp, sCode) = execCmdEx(
    BinPath & &" scatter -n 4 -o {tmpDir}/shard.vcf.gz {SmallVcf} 2>&1")
  doAssert sCode == 0, &"G06 scatter exited {sCode}:\n{sOutp}"
  # Collect shard paths.
  var shards: seq[string]
  for i in 1..4:
    shards.add(tmpDir / ("shard_" & $i & ".shard.vcf.gz"))
  let shardsArg = shards.join(" ")
  let outPath = tmpDir / "gathered.vcf.gz"
  let (gOutp, gCode) = runGatherSubcmd(&"-o {outPath} {shardsArg}")
  doAssert gCode == 0, &"G06 gather exited {gCode}:\n{gOutp}"
  doAssert fileExists(outPath), "G06: output missing"
  let got = countRecords(outPath)
  let orig = countRecords(SmallVcf)
  doAssert got == orig, &"G06: record count {got} != {orig}"
  doAssert recordsHash(outPath) == recordsHash(SmallVcf),
    "G06: content hash mismatch"
  removeDir(tmpDir)

# ---------------------------------------------------------------------------
# G07 — BCF: scatter → gather, record count and hash match
# ---------------------------------------------------------------------------

timed("G07", "gather subcommand BCF: records, content hash matches"):
  let tmpDir = createTempDir("blocky_", "")
  let (sOutp, sCode) = execCmdEx(
    BinPath & &" scatter -n 4 -o {tmpDir}/shard.bcf {SmallBcf} 2>&1")
  doAssert sCode == 0, &"G07 scatter exited {sCode}:\n{sOutp}"
  var shards: seq[string]
  for i in 1..4:
    shards.add(tmpDir / ("shard_" & $i & ".shard.bcf"))
  let shardsArg = shards.join(" ")
  let outPath = tmpDir / "gathered.bcf"
  let (gOutp, gCode) = runGatherSubcmd(&"-o {outPath} {shardsArg}")
  doAssert gCode == 0, &"G07 gather exited {gCode}:\n{gOutp}"
  doAssert fileExists(outPath), "G07: output missing"
  let got = countRecords(outPath)
  let orig = countRecords(SmallBcf)
  doAssert got == orig, &"G07: record count {got} != {orig}"
  doAssert recordsHash(outPath) == recordsHash(SmallBcf),
    "G07: content hash mismatch"
  removeDir(tmpDir)

# ---------------------------------------------------------------------------
# G08 — Omit -o: output goes to stdout (uncompressed VCF, record count matches)
# ---------------------------------------------------------------------------

timed("G08", "gather subcommand stdout: records written to stdout"):
  let tmpDir = createTempDir("blocky_", "")
  let (sOutp, sCode) = execCmdEx(
    BinPath & &" scatter -n 4 -o {tmpDir}/shard.vcf.gz {SmallVcf} 2>&1")
  doAssert sCode == 0, &"G08 scatter exited {sCode}:\n{sOutp}"
  var shards: seq[string]
  for i in 1..4:
    shards.add(tmpDir / ("shard_" & $i & ".shard.vcf.gz"))
  let shardsArg = shards.join(" ")
  # Capture stdout to a file via shell redirection; no -o flag.
  let stdoutFile = tmpDir / "stdout.vcf"
  let (gOutp, gCode) = execCmdEx(
    BinPath & &" gather {shardsArg} > {stdoutFile} 2>&1")
  doAssert gCode == 0, &"G08 gather stdout exited {gCode}:\n{gOutp}"
  doAssert fileExists(stdoutFile), "G08: stdout capture file missing"
  let got  = countRecords(stdoutFile)
  let orig = countRecords(SmallVcf)
  doAssert got == orig, &"G08: record count {got} != {orig}"
  removeDir(tmpDir)

# ---------------------------------------------------------------------------
# G09 — No input files exits non-zero
# ---------------------------------------------------------------------------

timed("G09", "gather subcommand: no input files exits non-zero"):
  let tmpDir = createTempDir("blocky_", "")
  let outPath = tmpDir / "out.vcf.gz"
  let (outp, code) = runGatherSubcmd(&"-o {outPath}")
  doAssert code != 0, "G09: no inputs should exit non-zero"
  removeDir(tmpDir)

# ---------------------------------------------------------------------------
# G10 — Missing input file exits non-zero
# ---------------------------------------------------------------------------

timed("G10", "gather subcommand: missing input file exits non-zero"):
  let tmpDir = createTempDir("blocky_", "")
  let outPath = tmpDir / "out.vcf.gz"
  let (outp, code) = runGatherSubcmd(&"-o {outPath} /nonexistent/shard.vcf.gz")
  doAssert code != 0, "G10: missing input file should exit non-zero"
  doAssert "not found" in outp.toLowerAscii,
    &"G10: error should mention 'not found', got: {outp}"
  removeDir(tmpDir)

# ===========================================================================
# G11–G17 — #CHROM header validation tests
# ===========================================================================

# ---------------------------------------------------------------------------
# G11 — extractChromLine: finds #CHROM in VCF header bytes
# ---------------------------------------------------------------------------

timed("G11", "extractChromLine: found #CHROM in VCF header bytes"):
  var header: seq[byte]
  for c in "##fileformat=VCFv4.2\n##INFO=<ID=DP>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n":
    header.add(byte(c))
  let line = extractChromLine(header)
  doAssert line.startsWith("#CHROM"), &"G11: expected '#CHROM...' got '{line}'"

# ---------------------------------------------------------------------------
# G12 — extractChromLine: returns "" when not found
# ---------------------------------------------------------------------------

timed("G12", "extractChromLine: returns empty string when not found"):
  var data: seq[byte]
  for c in "##fileformat=VCFv4.2\n##INFO=<ID=DP>\n":
    data.add(byte(c))
  let line = extractChromLine(data)
  doAssert line == "", &"G12: expected '' got '{line}'"

# ---------------------------------------------------------------------------
# G13 — chromLineFromBytes: works on BGZF VCF file bytes
# ---------------------------------------------------------------------------

timed("G13", "chromLineFromBytes: found #CHROM in BGZF VCF"):
  doAssert fileExists(SmallVcf), "G13: VCF fixture missing"
  let fileSize = getFileSize(SmallVcf).int
  var allBytes = newSeq[byte](fileSize)
  let f = open(SmallVcf, fmRead)
  discard readBytes(f, allBytes, 0, fileSize)
  f.close()
  let (fmt, isBgzf) = sniffStreamFormat(allBytes)
  doAssert fmt == ffVcf, &"G13: expected ffVcf, got {fmt}"
  let line = chromLineFromBytes(allBytes, fmt, isBgzf)
  doAssert line.startsWith("#CHROM"), &"G13: expected '#CHROM...' got '{line}'"

# ---------------------------------------------------------------------------
# G14 — chromLineFromBytes: works on BGZF BCF file bytes
# ---------------------------------------------------------------------------

timed("G14", "chromLineFromBytes: found #CHROM in BGZF BCF"):
  doAssert fileExists(SmallBcf), "G14: BCF fixture missing"
  let fileSize = getFileSize(SmallBcf).int
  var allBytes = newSeq[byte](fileSize)
  let f = open(SmallBcf, fmRead)
  discard readBytes(f, allBytes, 0, fileSize)
  f.close()
  let (fmt, isBgzf) = sniffStreamFormat(allBytes)
  doAssert fmt == ffBcf, &"G14: expected ffBcf, got {fmt}"
  let line = chromLineFromBytes(allBytes, fmt, isBgzf)
  doAssert line.startsWith("#CHROM"), &"G14: expected '#CHROM...' got '{line}'"

# ---------------------------------------------------------------------------
# G15 — chromLineFromFile: matches chromLineFromBytes
# ---------------------------------------------------------------------------

timed("G15", "chromLineFromFile: matches chromLineFromBytes result"):
  doAssert fileExists(SmallVcf), "G15: VCF fixture missing"
  let fileSize = getFileSize(SmallVcf).int
  var allBytes = newSeq[byte](fileSize)
  let f = open(SmallVcf, fmRead)
  discard readBytes(f, allBytes, 0, fileSize)
  f.close()
  let (fmt, isBgzf) = sniffStreamFormat(allBytes)
  let fromBytes = chromLineFromBytes(allBytes, fmt, isBgzf)
  let fromFile  = chromLineFromFile(SmallVcf, fmt, isBgzf)
  doAssert fromBytes == fromFile,
    &"G15: chromLineFromFile != chromLineFromBytes (got '{fromFile}' vs '{fromBytes}')"
  doAssert fromBytes.startsWith("#CHROM"), &"G15: expected '#CHROM...' got '{fromBytes}'"

# ---------------------------------------------------------------------------
# G16 — gather subcommand: matching #CHROM lines → success
# ---------------------------------------------------------------------------

timed("G16", "gather #CHROM match: success"):
  let tmpDir = createTempDir("blocky_", "")
  let (sOutp, sCode) = execCmdEx(
    BinPath & &" scatter -n 2 -o {tmpDir}/shard.vcf.gz {SmallVcf} 2>&1")
  doAssert sCode == 0, &"G16 scatter exited {sCode}:\n{sOutp}"
  var shards: seq[string]
  for i in 1..2:
    shards.add(tmpDir / ("shard_" & $i & ".shard.vcf.gz"))
  let outPath = tmpDir / "merged.vcf.gz"
  let (gOutp, gCode) = runGatherSubcmd(&"-o {outPath} " & shards.join(" "))
  doAssert gCode == 0, &"G16: gather with matching #CHROM should succeed, got {gCode}:\n{gOutp}"
  let got  = countRecords(outPath)
  let orig = countRecords(SmallVcf)
  doAssert got == orig, &"G16: record count {got} != {orig}"
  removeDir(tmpDir)

# ---------------------------------------------------------------------------
# G17 — gather subcommand: mismatched #CHROM → exit 1, no partial output
# ---------------------------------------------------------------------------

timed("G17", "gather #CHROM mismatch: exits 1, no partial output"):
  # Build two synthetic BGZF-VCF files with different sample columns.
  let tmpDir = createTempDir("blocky_", "")
  let vcfA =
    "##fileformat=VCFv4.2\n" &
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSampleA\n" &
    "1\t100\t.\tA\tT\t.\tPASS\t.\tGT\t0/1\n"
  let vcfB =
    "##fileformat=VCFv4.2\n" &
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSampleB\n" &
    "1\t200\t.\tG\tC\t.\tPASS\t.\tGT\t1/1\n"
  var rawA, rawB: seq[byte]
  for c in vcfA: rawA.add(byte(c))
  for c in vcfB: rawB.add(byte(c))
  let shardA = tmpDir / "shardA.vcf.gz"
  let shardB = tmpDir / "shardB.vcf.gz"
  block:
    var dataA = compressToBgzfMulti(rawA); dataA.add(@BGZF_EOF)
    let fA = open(shardA, fmWrite)
    discard fA.writeBytes(dataA, 0, dataA.len); fA.close()
  block:
    var dataB = compressToBgzfMulti(rawB); dataB.add(@BGZF_EOF)
    let fB = open(shardB, fmWrite)
    discard fB.writeBytes(dataB, 0, dataB.len); fB.close()

  let outPath = tmpDir / "merged.vcf.gz"
  let (gOutp, gCode) = runGatherSubcmd(&"-o {outPath} {shardA} {shardB}")
  doAssert gCode != 0, &"G17: mismatched #CHROM should exit non-zero, got {gCode}"
  doAssert "chrom" in gOutp.toLowerAscii or "mismatch" in gOutp.toLowerAscii,
    &"G17: error message should mention mismatch, got: {gOutp}"
  removeDir(tmpDir)
