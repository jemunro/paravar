## Tests for bgzf.nim — BGZF I/O (V1), format sniffing (V2), file format sniffing (V3).
## Run from project root: nim c -r tests/test_vcf_utils.nim

echo "--------------- Test BGZF ---------------"

import std/[algorithm, os, strformat]
import test_utils
import "../src/blocky/bgzf"

# libdeflate CRC32 — already linked via -ldeflate in bgzf
proc libdeflateCrc32(crc: cuint; buf: pointer; len: csize_t): cuint
  {.importc: "libdeflate_crc32", header: "<libdeflate.h>".}

proc dataCrc32(data: seq[byte]): uint32 =
  ## Compute CRC32 of data using libdeflate.
  if data.len == 0:
    return libdeflateCrc32(0'u32, nil, 0).uint32
  libdeflateCrc32(0'u32, data[0].unsafeAddr, data.len.csize_t).uint32

const DataDir  = "tests/data"
const SmallVcf = DataDir / "small.vcf.gz"
const TinyVcf  = DataDir / "tiny.vcf.gz"
const SmallBcf = DataDir / "small.bcf"

# ---------------------------------------------------------------------------
# Helper: read raw bytes from a file slice
# ---------------------------------------------------------------------------
proc leU32At(data: seq[byte]; pos: int): uint32 =
  data[pos].uint32 or (data[pos+1].uint32 shl 8) or
  (data[pos+2].uint32 shl 16) or (data[pos+3].uint32 shl 24)

proc readFileSlice(path: string; start: int64; length: int): seq[byte] =
  let f = open(path, fmRead)
  defer: f.close()
  f.setFilePos(start)
  result = newSeq[byte](length)
  discard readBytes(f, result, 0, length)

# ===========================================================================
# V1 — BGZF I/O: block scanning, raw copy, compress/decompress, split, CRC
# ===========================================================================

# ---------------------------------------------------------------------------
# V1.1 — scanBgzfBlockStarts: offsets valid, first=0, last block is 28-byte EOF
# ---------------------------------------------------------------------------

timed("V1.1", "scanBgzfBlockStarts: offsets valid, first=0, last is EOF"):
  let starts = scanBgzfBlockStarts(SmallVcf)
  doAssert starts.len >= 2,
    &"expected >= 2 blocks (data + EOF), got {starts.len}"
  doAssert starts[0] == 0,
    &"first block must start at offset 0, got {starts[0]}"
  for off in starts:
    let hdr = readFileSlice(SmallVcf, off, 3)
    doAssert hdr[0] == 0x1f and hdr[1] == 0x8b and hdr[2] == 0x08,
      &"bad BGZF magic at offset {off}"
  let lastHdr = readFileSlice(SmallVcf, starts[^1], 18)
  let blkSize = bgzfBlockSize(lastHdr)
  doAssert blkSize == 28,
    &"expected EOF block of 28 bytes, got {blkSize}"

# ---------------------------------------------------------------------------
# V1.2 — scanBgzfBlockStarts range: range-limited scan truncates correctly
# ---------------------------------------------------------------------------

timed("V1.2", "scanBgzfBlockStarts range: truncates at upper bound"):
  let allStarts = scanBgzfBlockStarts(SmallVcf)
  doAssert allStarts.len >= 2
  let first = scanBgzfBlockStarts(SmallVcf, 0, allStarts[1])
  doAssert first == @[allStarts[0]],
    &"range scan: expected [{allStarts[0]}], got {first}"

# ---------------------------------------------------------------------------
# V1.4 — compressToBgzf/decompressBgzf round-trip
# ---------------------------------------------------------------------------

timed("V1.4", "compressToBgzf/decompressBgzf: round-trip"):
  let original = "Hello, BGZF world!\nSecond line.\n"
  let origBytes = cast[seq[byte]](original)
  let compressed = compressToBgzf(origBytes)
  doAssert compressed[0] == 0x1f and compressed[1] == 0x8b,
    "compressed output missing gzip magic"
  doAssert compressed[12] == 0x42 and compressed[13] == 0x43,
    "compressed output missing BC extra field"
  let decompressed = decompressBgzf(compressed)
  doAssert decompressed == origBytes,
    &"round-trip mismatch: {decompressed} != {origBytes}"

# ---------------------------------------------------------------------------
# V1.5 — round-trip empty: empty input compresses and decompresses correctly
# ---------------------------------------------------------------------------

timed("V1.5", "round-trip empty: empty input"):
  let compressed = compressToBgzf(@[])
  doAssert bgzfBlockSize(compressed) > 0, "empty compress: invalid block header"
  let decompressed = decompressBgzf(compressed)
  doAssert decompressed.len == 0, "empty round-trip: expected empty result"

# ---------------------------------------------------------------------------
# V1.6 — decompressBgzf fixture: real fixture block starts with '#'
# ---------------------------------------------------------------------------

timed("V1.6", "decompressBgzf fixture: real block decompresses"):
  let starts = scanBgzfBlockStarts(SmallVcf)
  var dataStart = -1'i64
  var dataSize  = 0
  for i, off in starts:
    let hdr = readFileSlice(SmallVcf, off, 18)
    let sz = bgzfBlockSize(hdr)
    if sz != 28:
      dataStart = off
      dataSize  = sz
      break
  doAssert dataStart >= 0, "no non-EOF block found in small.vcf.gz"
  let raw = readFileSlice(SmallVcf, dataStart, dataSize)
  let decompressed = decompressBgzf(raw)
  doAssert decompressed.len > 0, "decompressed data block is empty"
  doAssert decompressed[0] == byte('#') or decompressed[0] == byte('c'),
    &"unexpected first byte: {decompressed[0]}"

# V1.7-V1.11 (splitChunk, bcfFirstDataOffset, splitBcfBoundaryBlock) removed:
# these procs no longer exist after Milestone V — scatter splits exclusively
# at index virtual offsets via splitBgzfBlockAtUOffset, eliminating the
# midpoint-line/record search and the bcfFirstDataOffset helper. End-to-end
# scatter correctness for both formats is covered by CL10-CL13 in test_cli.nim.

# ---------------------------------------------------------------------------
# V1.12 — BGZF CRC32 validation: stored CRC matches computed value
# ---------------------------------------------------------------------------

timed("V1.12", "BGZF CRC32: stored matches computed"):
  doAssert fileExists(SmallVcf), "fixture missing"
  let starts = scanBgzfBlockStarts(SmallVcf)
  var testOff = -1'i64
  for off in starts:
    let hdr = readFileSlice(SmallVcf, off, 18)
    let sz = bgzfBlockSize(hdr)
    if sz != 28:
      testOff = off
      break
  doAssert testOff >= 0, "no non-EOF block found in " & SmallVcf

  let hdr = readFileSlice(SmallVcf, testOff, 18)
  let blkSize = bgzfBlockSize(hdr)
  let blk = readFileSlice(SmallVcf, testOff, blkSize)

  let storedCrc = leU32At(blk, blkSize - 8)
  let decompressed = decompressBgzf(blk)
  let computedCrc = dataCrc32(decompressed)

  doAssert storedCrc != 0, "stored CRC32 should be non-zero for real data"
  doAssert storedCrc == computedCrc,
    &"BGZF CRC32 mismatch: stored={storedCrc:#x} computed={computedCrc:#x}"

# ===========================================================================
# V2 — Format sniffing: isBgzfStream, sniffFormat, sniffStreamFormat
# ===========================================================================

# ---------------------------------------------------------------------------
# V2.1 — isBgzfStream: BGZF magic detected; plain gzip and random bytes rejected
# ---------------------------------------------------------------------------

timed("V2.1", "isBgzfStream: BGZF magic detected, non-BGZF rejected"):
  let bgzfHead = [0x1f'u8, 0x8b, 0x08, 0x04, 0x00]
  doAssert isBgzfStream(bgzfHead), "should detect BGZF magic"

  let gzHead = [0x1f'u8, 0x8b, 0x08, 0x00, 0x00]
  doAssert not isBgzfStream(gzHead), "plain gzip is not BGZF"

  let rnd = [0x42'u8, 0x43, 0x46, 0x02]
  doAssert not isBgzfStream(rnd), "random bytes are not BGZF"

  let short = [0x1f'u8, 0x8b]
  doAssert not isBgzfStream(short), "too-short buffer is not BGZF"

# ---------------------------------------------------------------------------
# V2.2 — sniffFormat: BCF/VCF/text detected from uncompressed bytes
# ---------------------------------------------------------------------------

timed("V2.2", "sniffFormat: BCF/VCF/text detected correctly"):
  let bcfBytes = [byte('B'), byte('C'), byte('F'), 0x02'u8, 0x02'u8, 0x00'u8]
  doAssert sniffFormat(bcfBytes) == ffBcf, "BCF magic -> ffBcf"

  var vcfBytes: seq[byte]
  for c in "##fileformatVCFv4.2\n":
    vcfBytes.add(byte(c))
  doAssert sniffFormat(vcfBytes) == ffVcf, "##fileformat -> ffVcf"

  var txtBytes: seq[byte]
  for c in "CHROM\tPOS\tID\n":
    txtBytes.add(byte(c))
  doAssert sniffFormat(txtBytes) == ffText, "other bytes -> ffText"

  doAssert sniffFormat([0x00'u8, 0x01'u8]) == ffText, "short buffer -> ffText"

  let bcfExact = [byte('B'), byte('C'), byte('F'), 0x02'u8, 0x02'u8]
  doAssert sniffFormat(bcfExact) == ffBcf, "exact BCF magic length -> ffBcf"

# ---------------------------------------------------------------------------
# V2.3 — sniffStreamFormat BCF: small.bcf detected as BCF/BGZF
# ---------------------------------------------------------------------------

timed("V2.3", "sniffStreamFormat: BCF/BGZF"):
  doAssert fileExists(SmallBcf), "BCF fixture missing"
  let f = open(SmallBcf, fmRead)
  var buf = newSeq[byte](65536)
  let n = readBytes(f, buf, 0, 65536)
  f.close()
  buf.setLen(n)
  let (fmt, isBgzf) = sniffStreamFormat(buf)
  doAssert fmt == ffBcf,  &"small.bcf: expected ffBcf, got {fmt}"
  doAssert isBgzf,         "small.bcf: expected BGZF stream"

# ---------------------------------------------------------------------------
# V2.4 — sniffStreamFormat VCF: small.vcf.gz detected as VCF/BGZF
# ---------------------------------------------------------------------------

timed("V2.4", "sniffStreamFormat: VCF/BGZF"):
  doAssert fileExists(SmallVcf), "VCF fixture missing"
  let f = open(SmallVcf, fmRead)
  var buf = newSeq[byte](65536)
  let n = readBytes(f, buf, 0, 65536)
  f.close()
  buf.setLen(n)
  let (fmt, isBgzf) = sniffStreamFormat(buf)
  doAssert fmt == ffVcf,  &"small.vcf.gz: expected ffVcf, got {fmt}"
  doAssert isBgzf,         "small.vcf.gz: expected BGZF stream"

# ---------------------------------------------------------------------------
# V2.5 — sniffStreamFormat text: plain text detected as text/uncompressed
# ---------------------------------------------------------------------------

timed("V2.5", "sniffStreamFormat: text/uncompressed"):
  var raw: seq[byte]
  for c in "col1\tcol2\tcol3\nhello\tworld\t42\n":
    raw.add(byte(c))
  let (fmt, isBgzf) = sniffStreamFormat(raw)
  doAssert fmt == ffText, &"plain text: expected ffText, got {fmt}"
  doAssert not isBgzf,    "plain text: should not be BGZF"

# ---------------------------------------------------------------------------
# V2.6 — sniffStreamFormat uncompressed VCF: VCF/uncompressed
# ---------------------------------------------------------------------------

timed("V2.6", "sniffStreamFormat: uncompressed VCF"):
  var raw: seq[byte]
  for c in "##fileformatVCFv4.2\n##source=blocky\n#CHROM\tPOS\n":
    raw.add(byte(c))
  let (fmt, isBgzf) = sniffStreamFormat(raw)
  doAssert fmt == ffVcf,  &"uncompressed VCF: expected ffVcf, got {fmt}"
  doAssert not isBgzf,    "uncompressed VCF: should not be BGZF"

# ---------------------------------------------------------------------------
# V2.7 — sniffStreamFormat compressed text: BGZF text detected as text/BGZF
# ---------------------------------------------------------------------------

timed("V2.7", "sniffStreamFormat: BGZF-compressed text"):
  var raw: seq[byte]
  for c in "hello world\n":
    raw.add(byte(c))
  let compressed = compressToBgzf(raw)
  let (fmt, isBgzf) = sniffStreamFormat(compressed)
  doAssert fmt == ffText, &"BGZF text: expected ffText, got {fmt}"
  doAssert isBgzf,         "BGZF text: expected BGZF stream"

# ===========================================================================
# V3 — File format sniffing: sniffFileFormat
# ===========================================================================

# ---------------------------------------------------------------------------
# V3.1 — sniffFileFormat VCF BGZF: small.vcf.gz
# ---------------------------------------------------------------------------

timed("V3.1", "sniffFileFormat: VCF BGZF"):
  doAssert fileExists(SmallVcf), "fixture missing: " & SmallVcf
  let (fmt, compressed) = sniffFileFormat(SmallVcf)
  doAssert fmt == ffVcf,  "V3.1: expected ffVcf, got " & $fmt
  doAssert compressed,    "V3.1: expected compressed=true"

# ---------------------------------------------------------------------------
# V3.2 — sniffFileFormat BCF BGZF: small.bcf
# ---------------------------------------------------------------------------

timed("V3.2", "sniffFileFormat: BCF BGZF"):
  doAssert fileExists(SmallBcf), "fixture missing: " & SmallBcf
  let (fmt, compressed) = sniffFileFormat(SmallBcf)
  doAssert fmt == ffBcf,  "V3.2: expected ffBcf, got " & $fmt
  doAssert compressed,    "V3.2: expected compressed=true"

# ---------------------------------------------------------------------------
# V3.3 — sniffFileFormat uncompressed VCF
# ---------------------------------------------------------------------------

timed("V3.3", "sniffFileFormat: uncompressed VCF"):
  let tmp = getTempDir() / "blocky_sniff_test_v3_3.vcf"
  writeFile(tmp, "##fileformat=VCFv4.1\n#CHROM\tPOS\tID\tREF\tALT\n")
  let (fmt, compressed) = sniffFileFormat(tmp)
  doAssert fmt == ffVcf,   "V3.3: expected ffVcf, got " & $fmt
  doAssert not compressed, "V3.3: expected compressed=false"
  removeFile(tmp)

# ---------------------------------------------------------------------------
# V3.4 — sniffFileFormat uncompressed BCF (magic bytes only)
# ---------------------------------------------------------------------------

timed("V3.4", "sniffFileFormat: uncompressed BCF magic"):
  let tmp = getTempDir() / "blocky_sniff_test_v3_4.bcf_raw"
  let magic: seq[byte] = @[byte('B'), byte('C'), byte('F'), 0x02'u8, 0x02'u8]
  var f = open(tmp, fmWrite)
  discard f.writeBytes(magic, 0, magic.len)
  f.close()
  let (fmt, compressed) = sniffFileFormat(tmp)
  doAssert fmt == ffBcf,   "V3.4: expected ffBcf, got " & $fmt
  doAssert not compressed, "V3.4: expected compressed=false"
  removeFile(tmp)

# ===========================================================================
# V4 — GZI index parsing
# ===========================================================================

const Kg22Vcf = DataDir / "chr22_1kg.vcf.gz"
const Kg22Gzi = Kg22Vcf & ".gzi"

# ---------------------------------------------------------------------------
# V4.1 — parseGziBlockStarts: offsets non-empty, sorted, first is 0
# ---------------------------------------------------------------------------

timed("V4.1", "parseGziBlockStarts: offsets non-empty, sorted, first=0"):
  if fileExists(Kg22Gzi):
    let starts = parseGziBlockStarts(Kg22Gzi)
    doAssert starts.len > 0, "GZI should have at least one entry"
    doAssert starts[0] == 0, &"first GZI offset should be 0, got {starts[0]}"
    for i in 1 ..< starts.len:
      doAssert starts[i] > starts[i - 1],
        &"GZI offsets not sorted at index {i}: {starts[i-1]} >= {starts[i]}"
  else:
    echo "  [skipped — chr22_1kg.vcf.gz.gzi not found]"

# ---------------------------------------------------------------------------
# V4.2 — parseGziBlockStarts cross-validate: every GZI offset is a valid block
# ---------------------------------------------------------------------------

timed("V4.2", "parseGziBlockStarts: every offset is a valid BGZF block start"):
  if fileExists(Kg22Gzi) and fileExists(Kg22Vcf):
    let gziStarts = parseGziBlockStarts(Kg22Gzi)
    let scanStarts = scanBgzfBlockStarts(Kg22Vcf)
    # Every GZI offset must appear in scanBgzfBlockStarts results.
    var scanSet: seq[int64] = scanStarts
    scanSet.sort()
    for off in gziStarts:
      let idx = scanSet.binarySearch(off)
      doAssert idx >= 0,
        &"GZI offset {off} not found in scanBgzfBlockStarts results"
  else:
    echo "  [skipped — chr22_1kg fixtures not found]"

# ===========================================================================
# V5 — Streaming BGZF compress / decompress
# ===========================================================================

# ---------------------------------------------------------------------------
# V5.1 — bgzfCompressStream/bgzfDecompressStream: round-trip identity
# ---------------------------------------------------------------------------

timed("V5.1", "bgzfCompressStream/bgzfDecompressStream: round-trip"):
  let original = decompressBgzfFile(SmallVcf)
  doAssert original.len > 0, "fixture decompressed to empty"

  let tmpComp = getTempDir() / "blocky_test_compress_stream.gz"
  let tmpDecomp = getTempDir() / "blocky_test_decompress_stream.raw"
  let tmpRaw = getTempDir() / "blocky_test_raw_input.bin"

  # Write raw bytes to a temp file, then compress via stream.
  var fRaw = open(tmpRaw, fmWrite)
  discard fRaw.writeBytes(original, 0, original.len)
  fRaw.close()

  var fIn = open(tmpRaw, fmRead)
  var fOut = open(tmpComp, fmWrite)
  bgzfCompressStream(fIn, fOut)
  fOut.close()
  fIn.close()

  # Verify output is valid BGZF (check magic at start).
  let compData = readFile(tmpComp)
  doAssert compData.len > 28, "compressed output too small"
  doAssert compData[0].byte == 0x1f and compData[1].byte == 0x8b,
    "compressed output missing gzip magic"

  # Decompress back via stream.
  fIn = open(tmpComp, fmRead)
  fOut = open(tmpDecomp, fmWrite)
  bgzfDecompressStream(fIn, fOut)
  fOut.close()
  fIn.close()

  let roundTrip = cast[seq[byte]](readFile(tmpDecomp))
  doAssert roundTrip.len == original.len,
    &"round-trip size mismatch: {roundTrip.len} vs {original.len}"
  doAssert roundTrip == original, "round-trip content mismatch"

  removeFile(tmpRaw)
  removeFile(tmpComp)
  removeFile(tmpDecomp)

# ---------------------------------------------------------------------------
# V5.2 — bgzfCompressStream: empty input produces only EOF block
# ---------------------------------------------------------------------------

timed("V5.2", "bgzfCompressStream: empty input -> EOF block only"):
  let tmpIn = getTempDir() / "blocky_test_empty_in.bin"
  let tmpOut = getTempDir() / "blocky_test_empty_out.gz"
  writeFile(tmpIn, "")
  var fIn = open(tmpIn, fmRead)
  var fOut = open(tmpOut, fmWrite)
  bgzfCompressStream(fIn, fOut)
  fOut.close()
  fIn.close()

  let data = readFile(tmpOut)
  doAssert data.len == 28, &"empty compress: expected 28-byte EOF, got {data.len}"

  removeFile(tmpIn)
  removeFile(tmpOut)
