## scatter — split a bgzipped VCF into N roughly equal shards.
##
## Algorithm phases:
##   1. Read TBI/CSI index for coarse BGZF block offsets.
##   2. Extract VCF header via hts-nim; record first-data-block offset.
##   3. Optimise shard boundaries using weighted bisection + block validation.
##   4. Write shards: header + raw middle blocks + recompressed boundary blocks + EOF.

import std/[algorithm, atomics, cpuinfo, os, posix, sequtils, sets, strformat, strutils]
{.warning[Deprecated]: off.}
import std/threadpool
{.warning[Deprecated]: on.}
import vcf_utils

# ---------------------------------------------------------------------------
# Optional info-level logging (enabled by --info flag via main.nim)
# ---------------------------------------------------------------------------

var verbose* = false

template info(msg: string) =
  if verbose: stderr.writeLine "info: " & msg

# ---------------------------------------------------------------------------
# Output path helpers
# ---------------------------------------------------------------------------

proc shardOutputPath*(tmpl: string; shardIdx: int; nShards: int): string =
  ## Resolve the output path for shard shardIdx (0-based).
  ## If tmpl contains "{}": replace with zero-padded 1-based shard number.
  ## Otherwise: prepend "shard_NN." to the basename of tmpl, keep directory.
  let padded = align($(shardIdx + 1), len($nShards), '0')
  if "{}" in tmpl:
    return tmpl.replace("{}", padded)
  let dir  = tmpl.parentDir
  let base = tmpl.lastPathPart
  let prefixed = "shard_" & padded & "." & base
  result = if dir.len == 0: prefixed else: dir / prefixed

# FileFormat and Compression are imported from vcf_utils.

# ---------------------------------------------------------------------------
# Internal binary-reader helpers (little-endian reads that advance a cursor)
# ---------------------------------------------------------------------------

proc readLeU32(data: openArray[byte]; pos: var int): uint32 =
  ## Read a little-endian uint32 from data at pos, then advance pos by 4.
  result = data[pos].uint32 or (data[pos+1].uint32 shl 8) or
           (data[pos+2].uint32 shl 16) or (data[pos+3].uint32 shl 24)
  pos += 4

proc readLeI32(data: openArray[byte]; pos: var int): int32 =
  ## Read a little-endian int32 from data at pos, then advance pos by 4.
  cast[int32](readLeU32(data, pos))

proc readLeU64(data: openArray[byte]; pos: var int): uint64 =
  ## Read a little-endian uint64 from data at pos, then advance pos by 8.
  let lo = readLeU32(data, pos).uint64   # advances pos by 4
  let hi = readLeU32(data, pos).uint64   # advances pos by 4 again
  result = lo or (hi shl 32)

# ---------------------------------------------------------------------------
# Index parsing — TBI
# ---------------------------------------------------------------------------

proc parseTbiBlockStarts*(tbiPath: string): seq[int64] =
  ## Parse a .tbi tabix index and return sorted unique BGZF file offsets
  ## (chunk-begin virtual offsets >> 16) of all indexed blocks.
  let raw = decompressBgzfFile(tbiPath)
  if raw.len < 8 or raw[0] != byte('T') or raw[1] != byte('B') or
     raw[2] != byte('I') or raw[3] != 0x01:
    quit(&"parseTbiBlockStarts: not a valid TBI index: {tbiPath}", 1)
  var pos = 4
  let nRef = readLeI32(raw, pos).int  # pos → 8
  pos += 24                            # skip 6 int32s (format/cols/meta/skip)
  let lNm = readLeI32(raw, pos).int   # pos → 36
  pos += lNm                          # skip sequence-name block
  var starts: seq[int64]
  for _ in 0 ..< nRef:
    let nBin = readLeU32(raw, pos).int
    for _ in 0 ..< nBin:
      pos += 4  # skip bin_id (uint32)
      let nChunk = readLeU32(raw, pos).int
      for _ in 0 ..< nChunk:
        let beg    = readLeU64(raw, pos)
        let endOff = readLeU64(raw, pos)
        if endOff > beg:
          starts.add((beg shr 16).int64)
    let nIntv = readLeU32(raw, pos).int
    pos += 8 * nIntv  # skip linear index intervals
  starts.sort()
  result = starts.deduplicate(isSorted = true)

proc parseTbiVirtualOffsets*(tbiPath: string): seq[(int64, int)] =
  ## Parse a .tbi index and return sorted unique virtual offsets as
  ## (block_file_offset, within_block_offset) pairs for all chunk starts.
  let raw = decompressBgzfFile(tbiPath)
  if raw.len < 8 or raw[0] != byte('T') or raw[1] != byte('B') or
     raw[2] != byte('I') or raw[3] != 0x01:
    quit(&"parseTbiVirtualOffsets: not a valid TBI index: {tbiPath}", 1)
  var pos = 4
  let nRef = readLeI32(raw, pos).int
  pos += 24                            # skip 6 int32s
  let lNm = readLeI32(raw, pos).int
  pos += lNm                          # skip sequence-name block
  var offsets: seq[(int64, int)]
  for _ in 0 ..< nRef:
    let nBin = readLeU32(raw, pos).int
    for _ in 0 ..< nBin:
      pos += 4  # skip bin_id
      let nChunk = readLeU32(raw, pos).int
      for _ in 0 ..< nChunk:
        let beg    = readLeU64(raw, pos)
        let endOff = readLeU64(raw, pos)
        if endOff > beg:
          offsets.add(((beg shr 16).int64, (beg and 0xFFFF).int))
    let nIntv = readLeU32(raw, pos).int
    pos += 8 * nIntv  # skip linear index intervals
  offsets.sort(proc(a, b: (int64, int)): int =
    if a[0] != b[0]: cmp(a[0], b[0]) else: cmp(a[1], b[1]))
  var deduped: seq[(int64, int)]
  for i, v in offsets:
    if i == 0 or v != offsets[i - 1]:
      deduped.add(v)
  result = deduped

# ---------------------------------------------------------------------------
# Index parsing — CSI
# ---------------------------------------------------------------------------

proc parseCsiVirtualOffsets*(csiPath: string): seq[(int64, int)] =
  ## Parse a .csi index and return sorted unique virtual offsets as
  ## (block_file_offset, within_block_offset) pairs for all chunk starts.
  let raw = decompressBgzfFile(csiPath)
  if raw.len < 8 or raw[0] != byte('C') or raw[1] != byte('S') or
     raw[2] != byte('I') or raw[3] != 0x01:
    quit(&"parseCsiVirtualOffsets: not a valid CSI index: {csiPath}", 1)
  var pos = 4
  pos += 8               # skip min_shift (int32) + depth (int32)
  let lAux = readLeI32(raw, pos).int
  pos += lAux
  let nRef = readLeI32(raw, pos).int
  var offsets: seq[(int64, int)]
  for _ in 0 ..< nRef:
    let nBin = readLeI32(raw, pos).int
    for _ in 0 ..< nBin:
      pos += 4   # skip bin (uint32)
      pos += 8   # skip loffset (uint64)
      let nChunk = readLeI32(raw, pos).int
      for _ in 0 ..< nChunk:
        let beg    = readLeU64(raw, pos)
        let endOff = readLeU64(raw, pos)
        if endOff > beg:
          offsets.add(((beg shr 16).int64, (beg and 0xFFFF).int))
  offsets.sort(proc(a, b: (int64, int)): int =
    if a[0] != b[0]: cmp(a[0], b[0]) else: cmp(a[1], b[1]))
  var deduped: seq[(int64, int)]
  for i, v in offsets:
    if i == 0 or v != offsets[i - 1]:
      deduped.add(v)
  result = deduped

proc parseCsiBlockStarts*(csiPath: string): seq[int64] =
  ## Parse a .csi index and return sorted unique BGZF file offsets.
  let raw = decompressBgzfFile(csiPath)
  if raw.len < 8 or raw[0] != byte('C') or raw[1] != byte('S') or
     raw[2] != byte('I') or raw[3] != 0x01:
    quit(&"parseCsiBlockStarts: not a valid CSI index: {csiPath}", 1)
  var pos = 4
  pos += 8               # skip min_shift (int32) + depth (int32)
  let lAux = readLeI32(raw, pos).int
  pos += lAux            # skip aux data
  let nRef = readLeI32(raw, pos).int
  var starts: seq[int64]
  for _ in 0 ..< nRef:
    let nBin = readLeI32(raw, pos).int
    for _ in 0 ..< nBin:
      pos += 4   # skip bin (uint32)
      pos += 8   # skip loffset (uint64)
      let nChunk = readLeI32(raw, pos).int
      for _ in 0 ..< nChunk:
        let beg    = readLeU64(raw, pos)
        let endOff = readLeU64(raw, pos)
        if endOff > beg:
          starts.add((beg shr 16).int64)
  starts.sort()
  result = starts.deduplicate(isSorted = true)

# ---------------------------------------------------------------------------
# Public entry point — detect index type and parse
# ---------------------------------------------------------------------------

proc readIndexBlockStarts*(vcfPath: string): seq[int64] =
  ## Detect a .tbi or .csi index alongside vcfPath and return sorted unique
  ## BGZF file offsets. Aborts with an error message if no index is found.
  let tbi = vcfPath & ".tbi"
  let csi = vcfPath & ".csi"
  if fileExists(csi):
    result = parseCsiBlockStarts(csi)
    info(&"csi: read {result.len} block offsets from {csi}")
  elif fileExists(tbi):
    result = parseTbiBlockStarts(tbi)
    info(&"tbi: read {result.len} block offsets from {tbi}")
  else:
    stderr.writeLine &"error: no index found for {vcfPath} (expected {tbi} or {csi})"
    quit(1)

proc readIndexVirtualOffsets*(vcfPath: string): seq[(int64, int)] =
  ## Detect a .csi or .tbi index and return sorted unique virtual offsets as
  ## (block_file_offset, within_block_offset) pairs. Returns @[] if no index.
  let csi = vcfPath & ".csi"
  let tbi = vcfPath & ".tbi"
  if fileExists(csi):
    result = parseCsiVirtualOffsets(csi)
  elif fileExists(tbi):
    result = parseTbiVirtualOffsets(tbi)
  else:
    result = @[]

proc scanAllBlockStarts*(vcfPath: string; firstDataBlock: int64): seq[int64] =
  ## Scan the entire BGZF file and return file offsets of all data blocks
  ## (offset >= firstDataBlock, EOF block excluded).
  ## Used when no .tbi or .csi index is available.
  let fileSize  = getFileSize(vcfPath)
  let eofOffset = fileSize - 28
  let all = scanBgzfBlockStarts(vcfPath)
  result = all.filterIt(it >= firstDataBlock and it < eofOffset)
  info(&"scan: found {result.len} data blocks (scanned {all.len} total)")

# ---------------------------------------------------------------------------
# Phase 2 — extract header and first-data-block offset
# ---------------------------------------------------------------------------

proc extractBcfHeader*(path: string): seq[byte] =
  ## Extract the BCF header blob (5-byte magic + 4-byte l_text + l_text bytes)
  ## from path, recompressed as BGZF.  Verifies the BCF magic and calls quit(1)
  ## on any format error.  Uses compressToBgzfMulti to handle headers > 65536 bytes.
  let starts = scanBgzfBlockStarts(path)
  let f = open(path, fmRead)
  defer: f.close()
  var accum: seq[byte]
  var lText = -1'i64
  var headerSize = -1'i64   # 5 + 4 + l_text
  for off in starts:
    var hdr = newSeq[byte](18)
    f.setFilePos(off)
    if readBytes(f, hdr, 0, 18) < 18: break
    let blkSize = bgzfBlockSize(hdr)
    if blkSize <= 0: break
    var blk = newSeq[byte](blkSize)
    f.setFilePos(off)
    discard readBytes(f, blk, 0, blkSize)
    accum.add(decompressBgzf(blk))
    if lText < 0 and accum.len >= 9:
      for i in 0 ..< 5:
        if accum[i] != BCF_MAGIC[i]:
          quit(&"extractBcfHeader: {path}: invalid BCF magic", 1)
      var p = 5
      lText = readLeU32(accum, p).int64
      headerSize = 5'i64 + 4'i64 + lText
    if headerSize >= 0 and accum.len.int64 >= headerSize:
      return compressToBgzfMulti(accum[0 ..< headerSize.int])
  if headerSize < 0:
    quit(&"extractBcfHeader: {path}: file too short to read BCF header", 1)
  quit(&"extractBcfHeader: {path}: header claims l_text={lText} but file " &
       &"only provides {accum.len} bytes", 1)

proc blockHasData(content: openArray[byte]; prevEndedWithNewline: bool): bool =
  ## Return true if content contains at least one complete line not starting
  ## with '#'.  prevEndedWithNewline must be true if the previous BGZF block
  ## ended with '\n' (i.e. the first byte of content is the start of a new
  ## line); if false, the first partial line is a continuation from the
  ## previous block and is skipped.
  var lineStart = 0
  # Skip the partial first line when we are mid-line from the previous block.
  if not prevEndedWithNewline:
    var found = false
    for i in 0 ..< content.len:
      if content[i] == byte('\n'):
        lineStart = i + 1
        found = true
        break
    if not found:
      return false   # entire block is a line continuation — no complete lines
  # Check each complete line (terminated by '\n').
  for i in lineStart ..< content.len:
    if content[i] == byte('\n'):
      if i > lineStart and content[lineStart] != byte('#'):
        return true
      lineStart = i + 1
  # Do NOT check the partial last line: it may be the start of a header line
  # whose '#' character lives in the next block.
  return false

proc getHeaderAndFirstBlock*(vcfPath: string): (seq[byte], int64) =
  ## Scan BGZF blocks to collect all VCF header lines ('#' lines) and locate
  ## the first block containing data.  No htslib dependency — reads raw bytes.
  ## Returns (header recompressed as BGZF, first-data-block file offset).
  ## Handles long header lines spanning BGZF block boundaries correctly.
  let f = open(vcfPath, fmRead)
  defer: f.close()
  var hdrBuf = newSeq[byte](18)
  if readBytes(f, hdrBuf, 0, 18) < 18 or bgzfBlockSize(hdrBuf) < 0:
    stderr.writeLine &"error: no BGZF blocks found in {vcfPath}"
    quit(1)
  var headerBytes: seq[byte]
  var pos = 0'i64
  var prevEndedWithNewline = true   # start of file is always a line boundary
  while true:
    f.setFilePos(pos)
    if readBytes(f, hdrBuf, 0, 18) < 18: break
    let blkSize = bgzfBlockSize(hdrBuf)
    if blkSize <= 0: break
    var blk = newSeq[byte](blkSize)
    f.setFilePos(pos)
    discard readBytes(f, blk, 0, blkSize)
    let content = decompressBgzf(blk)
    if content.len > 0:
      if blockHasData(content, prevEndedWithNewline):
        # First data block — extract any leading '#' lines to complete the header.
        # If the previous block ended mid-line, the continuation bytes here are
        # part of a '#' line already started in headerBytes; include them up to '\n'.
        var i = 0
        if not prevEndedWithNewline:
          while i < content.len:
            headerBytes.add(content[i])
            if content[i] == byte('\n'):
              i += 1
              break
            i += 1
        var lineStart = i
        while i < content.len:
          if content[i] == byte('\n'):
            if i > lineStart and content[lineStart] == byte('#'):
              for j in lineStart .. i: headerBytes.add(content[j])
            lineStart = i + 1
          i += 1
        let compressedHeader = compressToBgzfMulti(headerBytes)
        info(&"header: {headerBytes.len} bytes uncompressed, " &
             &"{compressedHeader.len} bytes as BGZF")
        info(&"first data block at file offset {pos}")
        return (compressedHeader, pos)
      else:
        # Pure header block — collect all decompressed bytes.
        headerBytes.add(content)
        prevEndedWithNewline = content[^1] == byte('\n')
    pos += blkSize.int64
  # Edge case: file has no data records (header only).
  let compressedHeader = compressToBgzfMulti(headerBytes)
  result = (compressedHeader, pos)

# ---------------------------------------------------------------------------
# Phase 3 — shard boundary optimisation
# ---------------------------------------------------------------------------

proc getLengths*(starts: seq[int64]; fileSize: int64): seq[int64] =
  ## Compute the byte length of each BGZF block: starts[i+1] - starts[i],
  ## with the last block running to fileSize.
  result = newSeq[int64](starts.len)
  for i in 0 ..< starts.len - 1:
    result[i] = starts[i + 1] - starts[i]
  if starts.len > 0:
    result[^1] = fileSize - starts[^1]

proc partitionBoundaries*(lengths: seq[int64]; n: int): seq[int] =
  ## Return n-1 block indices that best partition lengths into n equal-size parts.
  ## Mirrors Python's bisect_left on cumulative sums.
  ## Calls quit(1) if there are fewer blocks than n-1.
  if lengths.len < n - 1:
    let need = n - 1
    stderr.writeLine &"error: only {lengths.len} BGZF blocks but need at least" &
      &" {need} to create {n} shards"
    quit(1)
  var cum = newSeq[int64](lengths.len)
  cum[0] = lengths[0]
  for i in 1 ..< lengths.len:
    cum[i] = cum[i - 1] + lengths[i]
  let total = cum[^1]
  result = newSeq[int](n - 1)
  for i in 0 ..< n - 1:
    let target = (int64(i + 1) * total) div n
    result[i] = lowerBound(cum, target)

proc isValidBoundary*(vcfPath: string; start: int64; length: int64): bool =
  ## Return true if the BGZF block at start contains at least two complete
  ## lines — i.e., there is a '\n' that is not the last byte.
  ## This mirrors Python's is_valid_boundary check.
  let f = open(vcfPath, fmRead)
  f.setFilePos(start)
  var raw = newSeq[byte](length)
  let nRead = readBytes(f, raw, 0, length.int)
  f.close()
  if nRead < length.int: return false
  let content = decompressBgzf(raw)
  for i in 0 ..< content.len - 1:
    if content[i] == byte('\n'):
      return true   # '\n' found that is not the last byte
  return false

proc isValidBcfBoundary(bcfPath: string; start: int64; length: int64): bool =
  ## Return true if the BGZF block at start contains at least two complete
  ## BCF records — i.e., it can be split at a record boundary.
  let (valid, _, _) = splitBcfBoundaryBlock(bcfPath, start, length)
  result = valid

proc optimiseBoundaries*(vcfPath: string; startsIn: seq[int64]; nShards: int;
                          nThreads: int = 1;
                          format: FileFormat = ffVcf): (seq[int], seq[int64], seq[int64]) =
  ## Find shard boundary block indices that produce roughly equal-size parts.
  ## Single-pass algorithm:
  ##   1. Compute initial partition from coarse block starts.
  ##   2. Scan each boundary block range for finer BGZF sub-blocks (one pass).
  ##   3. Repartition once with the enriched block list.
  ##   4. Validate; fix any invalid boundary by finding the nearest valid neighbour.
  let fileSize = getFileSize(vcfPath)
  var lengths  = getLengths(startsIn, fileSize)
  var bounds   = partitionBoundaries(lengths, nShards)
  info(&"optimise: {startsIn.len} coarse blocks, initial boundary indices: {bounds}")

  # Scan each boundary block range once for fine-grained sub-blocks.
  var allStarts = startsIn
  var seen = allStarts.toHashSet   # O(1) duplicate check replaces O(n) `notin seq`
  if nThreads > 1:
    var scanFVs: seq[FlowVar[seq[int64]]]
    for bi in bounds:
      let s = startsIn[bi]; let l = lengths[bi]
      scanFVs.add(spawn scanBgzfBlockStarts(vcfPath, s, s + l))
    for fv in scanFVs:
      for off in ^fv:
        if off notin seen:
          seen.incl(off)
          allStarts.add(off)
  else:
    for bi in bounds:
      let s = startsIn[bi]; let l = lengths[bi]
      for off in scanBgzfBlockStarts(vcfPath, s, s + l):
        if off notin seen:
          seen.incl(off)
          allStarts.add(off)

  allStarts.sort()
  allStarts = allStarts.deduplicate(isSorted = true)
  lengths = getLengths(allStarts, fileSize)
  bounds  = partitionBoundaries(lengths, nShards)
  info(&"optimise: {allStarts.len} fine blocks after scan, " &
       &"boundary offsets: {bounds.mapIt(allStarts[it])}")

  # Validate boundaries; fix any invalid one by scanning neighbours.
  for i in 0 ..< bounds.len:
    let isValid = proc(bi: int): bool =
      if format == ffBcf: isValidBcfBoundary(vcfPath, allStarts[bi], lengths[bi])
      else:             isValidBoundary(vcfPath, allStarts[bi], lengths[bi])
    if not isValid(bounds[i]):
      var fixed = false
      for delta in 1 ..< allStarts.len:
        for candidate in [bounds[i] + delta, bounds[i] - delta]:
          if candidate >= 0 and candidate < allStarts.len and isValid(candidate):
            bounds[i] = candidate
            fixed = true
            break
        if fixed: break
      if not fixed:
        stderr.writeLine "error: could not find valid boundary near shard " &
          $(i + 1) & "; try fewer shards"
        quit(1)

  info(&"optimise: done (single pass)")
  result = (bounds, allStarts, lengths)

# ---------------------------------------------------------------------------
# Phase 4 — shard writing helpers
# ---------------------------------------------------------------------------

type SplitPair = object
  ## Value-type wrapper for the (head, tail) pair returned by splitChunk.
  ## Using an object (not a tuple) keeps FlowVar transfer across threads clean.
  head: seq[byte]
  tail: seq[byte]

proc splitChunkPair(vcfPath: string; offset: int64; size: int64): SplitPair {.gcsafe.} =
  ## Thin gcsafe wrapper around splitChunk for use with spawn.
  let (h, t) = splitChunk(vcfPath, offset, size)
  SplitPair(head: h, tail: t)

proc splitBcfChunkPair(bcfPath: string; offset: int64; size: int64): SplitPair {.gcsafe.} =
  ## Thin gcsafe wrapper around splitBcfBoundaryBlock for use with spawn.
  let (_, h, t) = splitBcfBoundaryBlock(bcfPath, offset, size)
  SplitPair(head: h, tail: t)

type ShardTask* = object
  ## All data needed to write one output shard.  Passed by value to worker
  ## threads so each thread owns its own copy of the small prepend/boundary
  ## buffers; the bulk of the data is read directly from vcfPath.
  vcfPath:      string
  outFd*:       cint        ## open writable file descriptor (file or pipe write-end)
  prepend:      seq[byte]   ## recompressed header + first-block tail / boundary tail
  rawStart:     int64       ## start of middle raw-copy region
  rawEnd:       int64       ## end of middle raw-copy region (exclusive)
  boundaryHead: seq[byte]   ## head half of boundary split; empty for last shard
  eofSeq:       seq[byte]   ## BGZF EOF block (same 28 bytes for every shard)
  logLine:      string      ## pre-formatted info line; printed when writing begins
  decompress*:  bool        ## if true, decompress BGZF blocks before writing

proc doWriteShard*(task: ShardTask): int {.gcsafe.} =
  ## Write one shard to its output fd.  Returns 0.  Designed for use with
  ## spawn so that multiple shards can be written in parallel.
  ## IOError (broken pipe) is caught silently: when writing to a pipe whose
  ## read-end has been closed by an early-exiting child, the failure is
  ## detected via waitpid's exit code rather than via a write error.
  ## When task.decompress is true, BGZF blocks are decompressed before writing
  ## so the receiver gets a raw uncompressed byte stream.
  if task.logLine.len > 0:
    stderr.writeLine "info: " & task.logLine
  var f: File
  discard open(f, FileHandle(task.outFd), fmWrite)
  try:
    if task.decompress:
      let prependRaw = decompressBgzfBytes(task.prepend)
      discard f.writeBytes(prependRaw, 0, prependRaw.len)
      if task.rawEnd > task.rawStart:
        decompressCopyBytes(task.vcfPath, f, task.rawStart,
                            task.rawEnd - task.rawStart)
      if task.boundaryHead.len > 0:
        let headRaw = decompressBgzfBytes(task.boundaryHead)
        discard f.writeBytes(headRaw, 0, headRaw.len)
      # No BGZF EOF block — the receiver gets a plain byte stream.
    else:
      discard f.writeBytes(task.prepend, 0, task.prepend.len)
      if task.rawEnd > task.rawStart:
        rawCopyBytes(task.vcfPath, f, task.rawStart, task.rawEnd - task.rawStart)
      if task.boundaryHead.len > 0:
        discard f.writeBytes(task.boundaryHead, 0, task.boundaryHead.len)
      discard f.writeBytes(task.eofSeq, 0, task.eofSeq.len)
  except IOError:
    discard  # broken pipe — child exited early; failure detected via waitpid
  try: f.close() except IOError: discard
  return 0

# ---------------------------------------------------------------------------
# Phase 4 — compute shard data and scatter entry point
# ---------------------------------------------------------------------------

proc computeShards*(vcfPath: string; nShards: int; nThreads: int = 1;
                    forceScan: bool = false;
                    format: FileFormat = ffVcf): seq[ShardTask] =
  ## Compute per-shard byte ranges and prepend buffers for vcfPath.
  ## Returns nShards ShardTask objects with outFd = -1 and logLine = "".
  ## Caller must set task.outFd before calling doWriteShard.
  ## nThreads controls threadpool parallelism for scanning and boundary splits.
  let fileSize = getFileSize(vcfPath)

  var headerBytes: seq[byte]
  var firstBlock: int64
  var starts: seq[int64]

  if format == ffBcf:
    # BCF: CSI index required; no TBI fallback; no auto-scan.
    # Use full CSI virtual offsets (block_off, u_off) for record-exact splitting.
    let csi = vcfPath & ".csi"
    if not fileExists(csi):
      stderr.writeLine &"error: BCF input requires a CSI index: {csi} not found"
      stderr.writeLine &"  (create one with 'bcftools index {vcfPath}')"
      quit(1)
    headerBytes = extractBcfHeader(vcfPath)
    let (firstBlockOff, firstUOff) = bcfFirstDataVirtualOffset(vcfPath)
    var voffs = parseCsiVirtualOffsets(csi)
    let firstVO = (firstBlockOff, firstUOff)
    if firstVO notin voffs: voffs.add(firstVO)
    voffs.sort(proc(a, b: (int64, int)): int =
      if a[0] != b[0]: cmp(a[0], b[0]) else: cmp(a[1], b[1]))
    var deduped: seq[(int64, int)]
    for i, v in voffs:
      if i == 0 or v != voffs[i - 1]: deduped.add(v)
    voffs = deduped
    info(&"bcf: {voffs.len} virtual offsets, firstData=({firstBlockOff},{firstUOff})")

    # Select nShards-1 boundary virtual offsets, evenly spaced by index.
    var boundaryVoffs: seq[(int64, int)]
    for i in 1 ..< nShards:
      let idx = (i * voffs.len) div nShards
      boundaryVoffs.add(voffs[idx])

    let eofStart = fileSize - 28
    let eofSeq: seq[byte] = @BGZF_EOF

    result = newSeq[ShardTask](nShards)
    for i in 0 ..< nShards:
      var prepend: seq[byte]
      prepend.add(headerBytes)

      # Determine the start virtual offset for this shard.
      let (startBlockOff, startUOff) =
        if i == 0: (firstBlockOff, firstUOff) else: boundaryVoffs[i - 1]

      # If start block is mid-record (startUOff > 0), prepend the tail of
      # that block (bytes [startUOff, end]).  If startUOff == 0, the block
      # starts cleanly and is included in the raw-copy region.
      if startUOff > 0:
        let (_, tail) = splitBgzfBlockAtUOffset(vcfPath, startBlockOff, startUOff)
        prepend.add(tail)

      # rawStart: if startUOff > 0, skip past the split block; if 0, include it.
      var hdrBuf = newSeq[byte](18)
      block getSize:
        let fTmp = open(vcfPath, fmRead)
        fTmp.setFilePos(startBlockOff)
        discard readBytes(fTmp, hdrBuf, 0, 18)
        fTmp.close()
      let startBlkSize = bgzfBlockSize(hdrBuf).int64
      let rawStart: int64 =
        if startUOff > 0: startBlockOff + startBlkSize else: startBlockOff

      var rawEnd: int64
      var boundaryHead: seq[byte]
      if i < nShards - 1:
        let (endBlockOff, endUOff) = boundaryVoffs[i]
        rawEnd = endBlockOff
        if endUOff > 0:
          let (head, _) = splitBgzfBlockAtUOffset(vcfPath, endBlockOff, endUOff)
          boundaryHead = head
      else:
        rawEnd = eofStart

      if rawEnd < rawStart: rawEnd = rawStart
      result[i] = ShardTask(vcfPath: vcfPath, outFd: -1, prepend: prepend,
                            rawStart: rawStart, rawEnd: rawEnd,
                            boundaryHead: boundaryHead, eofSeq: eofSeq,
                            logLine: "")
    return result  # BCF path complete; skip VCF path below

  else:
    # VCF: TBI or CSI index, or auto-scan fallback.
    let (hb, fb) = getHeaderAndFirstBlock(vcfPath)
    headerBytes = hb
    firstBlock  = fb
    let tbi      = vcfPath & ".tbi"
    let csi      = vcfPath & ".csi"
    let hasIndex = fileExists(csi) or fileExists(tbi)
    if forceScan or not hasIndex:
      if not hasIndex:
        stderr.writeLine &"warning: no index found for {vcfPath} — scanning BGZF blocks directly"
        stderr.writeLine "  (create an index with 'tabix -p vcf' for faster operation)"
      else:
        info("--force-scan: ignoring index, scanning all BGZF blocks")
      starts = scanAllBlockStarts(vcfPath, firstBlock)
    else:
      starts = readIndexBlockStarts(vcfPath)
      if firstBlock notin starts: starts.add(firstBlock)
      starts.sort()
      if starts.len >= 2:
        for off in scanBgzfBlockStarts(vcfPath, starts[0], starts[1]):
          if off notin starts: starts.add(off)
        starts.sort()

  let (bounds, finalStarts, lengths) =
    optimiseBoundaries(vcfPath, starts, nShards, nThreads, format = format)

  var splits: seq[SplitPair]
  if format == ffBcf:
    if nThreads > 1 and bounds.len > 0:
      var splitFVs = newSeq[FlowVar[SplitPair]](bounds.len)
      for idx, bi in bounds:
        splitFVs[idx] = spawn splitBcfChunkPair(vcfPath, finalStarts[bi], lengths[bi])
      for fv in splitFVs:
        splits.add(^fv)
    else:
      for bi in bounds:
        splits.add(splitBcfChunkPair(vcfPath, finalStarts[bi], lengths[bi]))
  else:
    if nThreads > 1 and bounds.len > 0:
      var splitFVs = newSeq[FlowVar[SplitPair]](bounds.len)
      for idx, bi in bounds:
        splitFVs[idx] = spawn splitChunkPair(vcfPath, finalStarts[bi], lengths[bi])
      for fv in splitFVs:
        splits.add(^fv)
    else:
      for bi in bounds:
        let (h, t) = splitChunk(vcfPath, finalStarts[bi], lengths[bi])
        splits.add(SplitPair(head: h, tail: t))

  let eofStart = fileSize - 28
  let eofSeq: seq[byte] = @BGZF_EOF
  info(&"computeShards: {nShards} shards, file size {fileSize} bytes, EOF at {eofStart}")

  result = newSeq[ShardTask](nShards)
  for i in 0 ..< nShards:
    var prepend: seq[byte]
    prepend.add(headerBytes)
    if i == 0:
      if (bounds.len == 0 or finalStarts[0] < finalStarts[bounds[0]]) and
         finalStarts.len >= 2:
        if format == ffBcf:
          prepend.add(removeBcfHeaderBytes(vcfPath))
        else:
          prepend.add(removeHeaderLines(vcfPath, 0, finalStarts[1]))
    else:
      prepend.add(splits[i - 1].tail)
    let rawStart: int64 =
      if i == 0:
        if finalStarts.len >= 2: finalStarts[1] else: eofStart
      else:
        finalStarts[bounds[i - 1] + 1]
    var rawEnd: int64 =
      if i == nShards - 1: eofStart else: finalStarts[bounds[i]]
    if rawEnd < rawStart: rawEnd = rawStart
    let boundaryHead = if i < nShards - 1: splits[i].head else: @[]
    result[i] = ShardTask(vcfPath: vcfPath, outFd: -1, prepend: prepend,
                          rawStart: rawStart, rawEnd: rawEnd,
                          boundaryHead: boundaryHead, eofSeq: eofSeq,
                          logLine: "")

# ---------------------------------------------------------------------------
# Interleaved block assignment
# ---------------------------------------------------------------------------

proc interleavedBlockAssignment*(nBlocks, nShards, chunkSize: int): seq[seq[Slice[int]]] =
  ## Assign nBlocks data blocks to nShards shards in round-robin chunks of
  ## chunkSize.  Returns a seq of length nShards; each element is a list of
  ## Slice[int] index ranges into the block array.  Every block index in
  ## 0..<nBlocks appears in exactly one shard's list.
  result = newSeq[seq[Slice[int]]](nShards)
  var pos = 0
  var chunkIdx = 0
  while pos < nBlocks:
    let shard = chunkIdx mod nShards
    let hi = min(pos + chunkSize, nBlocks) - 1  # inclusive end
    result[shard].add(pos .. hi)
    pos = hi + 1
    inc chunkIdx

# ---------------------------------------------------------------------------
# Interleaved shard writer
# ---------------------------------------------------------------------------

type InboxSlot* = object
  ## Single-slot mailbox for cross-worker partial-record handoff.
  ## Writer deposits buf+len then sets ready=true; reader spins on ready.
  ## ready uses atomic operations for cross-thread visibility.
  buf*:   array[4 * 1024 * 1024, byte]
  len*:   int32
  ready*: Atomic[bool]

type InboxArray* = object
  ## Array of per-worker inbox slots, heap-allocated for GC safety.
  slots*: ptr UncheckedArray[InboxSlot]
  count*: int

proc newInboxArray*(n: int): InboxArray =
  ## Allocate n inbox slots, all zeroed (ready = false).
  let p = cast[ptr UncheckedArray[InboxSlot]](allocShared0(n * sizeof(InboxSlot)))
  InboxArray(slots: p, count: n)

proc freeInboxArray*(inboxes: InboxArray) =
  if inboxes.slots != nil:
    deallocShared(inboxes.slots)

proc deposit*(inboxes: InboxArray; workerIdx: int; data: openArray[byte]) =
  ## Deposit data into workerIdx's inbox. Spins until the slot is free (ready ==
  ## false), then writes data and sets ready = true as a release signal.
  let slot = addr inboxes.slots[workerIdx]
  # Backpressure: wait for previous drain to complete before overwriting.
  while slot.ready.load(moAcquire):
    discard  # spin — slot still occupied by previous cycle
  let n = min(data.len, slot.buf.len)
  if n > 0:
    copyMem(addr slot.buf[0], unsafeAddr data[0], n)
  slot.len = n.int32
  slot.ready.store(true, moRelease)

proc drain*(inboxes: InboxArray; workerIdx: int; timeoutMs: int = 10000): seq[byte] =
  ## Spin-wait until workerIdx's inbox is ready, then return the data and reset.
  ## Raises an error if timeout is reached (depositing worker likely died).
  let slot = addr inboxes.slots[workerIdx]
  var waited = 0
  while not slot.ready.load(moAcquire):
    sleep(1)
    waited += 1
    if waited >= timeoutMs:
      stderr.writeLine "error: inbox drain timeout for worker " & $workerIdx &
        " — depositing worker may have crashed (waited " & $waited & "ms)"
      quit(1)
  result = newSeq[byte](slot.len)
  if slot.len > 0:
    copyMem(addr result[0], addr slot.buf[0], slot.len)
  slot.len = 0
  slot.ready.store(false, moRelease)

proc depositEmpty*(inboxes: InboxArray; workerIdx: int) =
  ## Signal that no handoff is needed (record boundary fell cleanly).
  ## Spins until the slot is free before signalling.
  let slot = addr inboxes.slots[workerIdx]
  while slot.ready.load(moAcquire):
    discard  # spin — slot still occupied by previous cycle
  slot.len = 0
  slot.ready.store(true, moRelease)

type InterleavedTask* = object
  vcfPath*:      string
  outFd*:        cint
  headerBytes*:  seq[byte]        ## uncompressed VCF/BCF header
  blockStarts*:  ptr seq[int64]   ## shared: all data block file offsets
  blockSizes*:   ptr seq[int64]   ## shared: all data block byte sizes
  chunkIndices*: seq[Slice[int]]  ## this shard's index ranges into blockStarts
  format*:       FileFormat       ## ffVcf or ffBcf
  csiVoffs*:     ptr seq[(int64, int)]  ## shared: CSI virtual offsets (BCF only; nil for VCF)
  shardIdx*:     int              ## this worker's index (0-based)
  nShards*:      int              ## total number of workers
  chunkSize*:    int              ## blocks per chunk (for chunkOwner calculation)
  inboxes*:      ptr InboxArray   ## shared inbox array

proc readAndDecompressBlocks(vcfPath: string; starts: ptr seq[int64];
                              sizes: ptr seq[int64];
                              slice: Slice[int]): seq[byte] {.gcsafe.} =
  ## Read and decompress all BGZF blocks in starts[slice] from vcfPath.
  ## Each starts[i]..starts[i]+sizes[i] range may contain multiple BGZF blocks
  ## (when index entries span multiple blocks); decompress each one in sequence.
  let f = open(vcfPath, fmRead)
  defer: f.close()
  result = @[]
  for i in slice:
    let offset = starts[][i]
    let sz     = sizes[][i].int
    var raw = newSeq[byte](sz)
    f.setFilePos(offset)
    discard readBytes(f, raw, 0, sz)
    var pos = 0
    while pos < raw.len:
      let blkSize = bgzfBlockSize(raw.toOpenArray(pos, raw.high))
      if blkSize <= 0: break
      result.add(decompressBgzf(raw.toOpenArray(pos, pos + blkSize - 1)))
      pos += blkSize

proc readAndDecompressOneBlock(vcfPath: string; start: int64;
                                size: int64): seq[byte] {.gcsafe.} =
  ## Read and decompress a single BGZF block.
  let f = open(vcfPath, fmRead)
  defer: f.close()
  var blk = newSeq[byte](size.int)
  f.setFilePos(start)
  discard readBytes(f, blk, 0, size.int)
  decompressBgzf(blk)

proc findFirstNewline(data: openArray[byte]): int =
  ## Return index of first '\n' in data, or -1 if none.
  for i in 0 ..< data.len:
    if data[i] == byte('\n'): return i
  return -1

proc findLastNewline(data: openArray[byte]): int =
  ## Return index of last '\n' in data, or -1 if none.
  for i in countdown(data.len - 1, 0):
    if data[i] == byte('\n'): return i
  return -1

proc findBcfRecordEnd*(data: openArray[byte]): int =
  ## Walk BCF records from the start of data.  Return the index just past the
  ## last complete record, or 0 if not even one complete record fits.
  var pos = 0
  while pos + 8 <= data.len:
    let lShared = cast[ptr uint32](unsafeAddr data[pos])[]
    let lIndiv  = cast[ptr uint32](unsafeAddr data[pos + 4])[]
    let recLen  = 8 + lShared.int + lIndiv.int
    if pos + recLen > data.len:
      break  # incomplete record
    pos += recLen
  return pos

proc headSkipFromIndex*(blockOff: int64; voffs: ptr seq[(int64, int)]): int =
  ## Look up the index virtual offsets to find the head skip for a block.
  ## Returns the u_off of the virtual offset entry with matching block_off.
  ## Used for both BCF (CSI) and indexed VCF (TBI/CSI). The caller guarantees
  ## blockOff was placed at an index entry boundary, so the lookup always
  ## succeeds when the index is well-formed.
  let v = voffs[]
  var lo = 0
  var hi = v.len - 1
  while lo <= hi:
    let mid = (lo + hi) div 2
    if v[mid][0] < blockOff:
      lo = mid + 1
    elif v[mid][0] > blockOff:
      hi = mid - 1
    else:
      return v[mid][1]
  return 0  # not found — caller's blockOff wasn't at an index boundary

proc chunkOwner*(blockIdx, nShards, chunkSize: int): int {.inline.} =
  ## Return the worker index that owns the chunk containing blockIdx.
  (blockIdx div chunkSize) mod nShards

proc writeInterleavedShard*(task: InterleavedTask): int {.gcsafe.} =
  ## Write one interleaved shard using the inbox model.
  ## For each chunk: deposit head to previous worker's inbox, write complete
  ## records, then unconditionally drain own inbox.  Returns 0.
  let nTotalBlocks = task.blockStarts[].len
  var f: File
  discard open(f, FileHandle(task.outFd), fmWrite)
  var completedChunks = 0
  try:
    # Write header
    if task.headerBytes.len > 0:
      discard f.writeBytes(task.headerBytes, 0, task.headerBytes.len)

    for ci, slice in task.chunkIndices:
      let buf = readAndDecompressBlocks(task.vcfPath, task.blockStarts,
                                        task.blockSizes, slice)

      # --- Head split: bytes before first record boundary ---
      var headEnd = 0
      if slice.a > 0:
        if task.csiVoffs != nil:
          # Indexed (BCF or VCF with TBI/CSI): exact u_off lookup
          headEnd = headSkipFromIndex(task.blockStarts[][slice.a], task.csiVoffs)
        else:
          # Unindexed VCF fallback: scan for first newline
          let nlPos = findFirstNewline(buf)
          if nlPos >= 0:
            headEnd = nlPos + 1
          else:
            headEnd = buf.len  # entire chunk is one partial line

        let prevOwner = chunkOwner(slice.a - 1, task.nShards, task.chunkSize)
        if prevOwner == task.shardIdx:
          # Previous chunk is ours — write head directly (no inbox needed).
          if headEnd > 0:
            discard f.writeBytes(buf, 0, headEnd)
        else:
          # Deposit to the other worker's inbox.
          if headEnd > 0:
            task.inboxes[].deposit(prevOwner, buf.toOpenArray(0, headEnd - 1))
          else:
            task.inboxes[].depositEmpty(prevOwner)

      # --- Write records from buf[headEnd..] ---
      # Write complete records, then trailing partial + inbox data.
      if headEnd < buf.len:
        let records = buf[headEnd ..< buf.len]
        var recEnd: int
        if task.format == ffBcf:
          recEnd = findBcfRecordEnd(records)
        else:
          let nlPos = findLastNewline(records)
          recEnd = if nlPos >= 0: nlPos + 1 else: 0
        if recEnd > 0:
          discard f.writeBytes(records, 0, recEnd)
        # Write trailing partial (bytes after last complete record).
        if recEnd < records.len:
          discard f.writeBytes(records, recEnd, records.len - recEnd)

      # --- Drain own inbox if the next chunk belongs to a different worker ---
      if slice.b + 1 < nTotalBlocks:
        let nextOwner = chunkOwner(slice.b + 1, task.nShards, task.chunkSize)
        if nextOwner != task.shardIdx:
          let inboxData = task.inboxes[].drain(task.shardIdx)
          if inboxData.len > 0:
            discard f.writeBytes(inboxData, 0, inboxData.len)

      completedChunks = ci + 1

  except IOError:
    discard  # broken pipe — child exited early

  # Fulfill inbox obligations for any remaining unprocessed chunks so other
  # workers don't deadlock waiting for deposits that will never come.
  for ci in completedChunks ..< task.chunkIndices.len:
    let slice = task.chunkIndices[ci]
    if slice.a > 0:
      let prevOwner = chunkOwner(slice.a - 1, task.nShards, task.chunkSize)
      if prevOwner != task.shardIdx:
        task.inboxes[].depositEmpty(prevOwner)
    # Also drain our own inbox if expected, so the depositing worker isn't stuck.
    if slice.b + 1 < nTotalBlocks:
      let nextOwner = chunkOwner(slice.b + 1, task.nShards, task.chunkSize)
      if nextOwner != task.shardIdx:
        discard task.inboxes[].drain(task.shardIdx)

  try: f.close() except IOError: discard
  return 0

proc scatter*(vcfPath: string; nShards: int; outputTemplate: string;
              nThreads: int = 1; forceScan: bool = false;
              format: FileFormat = ffVcf) =
  ## Split vcfPath into nShards bgzipped files.
  ## outputTemplate may contain {} (replaced with zero-padded shard number)
  ## or not (shard_NN. is prepended to the basename). mkdir -p is applied.
  ## nThreads controls parallelism; pass 0 to use all CPUs.
  ## Output is always BGZF compressed.
  let actualThreads = if nThreads == 0: countProcessors() else: nThreads
  setMaxPoolSize(actualThreads)
  info(&"scatter: using {actualThreads} thread(s)")
  var tasks = computeShards(vcfPath, nShards, actualThreads, forceScan, format)
  for i in 0 ..< nShards:
    let outPath = shardOutputPath(outputTemplate, i, nShards)
    createDir(outPath.parentDir)
    tasks[i].outFd = posix.open(outPath.cstring,
                                O_WRONLY or O_CREAT or O_TRUNC,
                                0o666.Mode)
    if tasks[i].outFd < 0:
      stderr.writeLine &"error: could not create output file: {outPath}"
      quit(1)
    if verbose:
      let rs = tasks[i].rawStart; let re = tasks[i].rawEnd
      let bn = if i < nShards - 1: &", boundary head {tasks[i].boundaryHead.len} bytes"
               else: ""
      tasks[i].logLine = &"shard {i+1}/{nShards}: {outPath} " &
        &"(prepend {tasks[i].prepend.len}B, raw {rs}..{re} = {re-rs}B{bn})"
  if actualThreads > 1:
    var writeFVs = newSeq[FlowVar[int]](nShards)
    for i, task in tasks:
      writeFVs[i] = spawn doWriteShard(task)
    for fv in writeFVs:
      discard ^fv
  else:
    for task in tasks:
      discard doWriteShard(task)
