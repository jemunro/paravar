# CLAUDE.md — paravar

Agentic development instructions for Claude Code. Read this file in full before starting any task.

---

## Project overview

**paravar** is a Nim CLI tool for parallelising VCF/BCF processing by exploiting BGZF block boundaries and tabix/CSI indexes.

The canonical reference for scatter behaviour is `example/scatter_vcf.py`.

---

## Subcommands

| Subcommand | Status | Description |
|---|---|---|
| `scatter` | ✅ | Split VCF/BCF into N shards |
| `run` | ✅ | Scatter + parallel tool pipelines → per-shard outputs |
| `run --gather` | ✅ | As above, gathered into single output file |
| `gather` | ✅ | Concatenate existing shard files into single output |

---

## CLI specification

### `scatter`

```
paravar scatter -n <n> -o <o> [options] <input>
```

Split input into N shards. Each shard is a valid standalone file.

### `run`

```
paravar run -n <n> -o <o> [options] <input> (--- | :::) <cmd> [(--- | :::) <cmd> ...]
```

Scatter + pipe each shard through a tool pipeline. Without `--gather`, writes N shard output files.

### `run --gather`

```
paravar run --gather [-o <o>] -n <n> [options] <input> (--- | :::) <cmd> [(--- | :::) <cmd> ...]
```

As `run`, but intercept each shard's stdout, strip headers from shards 2..N, and concatenate into `-o`. If `-o` is omitted or set to `/dev/stdout`, output is written to stdout.

### `gather`

```
paravar gather [-o <o>] [options] <shard1> [<shard2> ...]
```

Concatenate pre-existing shard files (output of `scatter` or `run`) into a single output file. Same header stripping and recompression logic as `run --gather`. No `--tmp-dir` needed — operates directly on input files. If `-o` is omitted or set to `/dev/stdout`, output is written to stdout.

---

## Pipeline separator

`---` (three dashes) and `:::` are both valid pipeline stage separators and may be mixed freely in the same invocation. `:::` is provided for familiarity with GNU parallel. Both split the argv into stages that are joined with `|` and executed via `sh -c`.

```bash
# These are equivalent:
paravar run -n 8 -o out.vcf.gz input.vcf.gz \
  --- bcftools view -i "GT='alt'" -Ou \
  --- bcftools view -s Sample -Oz

paravar run -n 8 -o out.vcf.gz input.vcf.gz \
  ::: bcftools view -i "GT='alt'" -Ou \
  ::: bcftools view -s Sample -Oz

# Mixing is allowed:
paravar run -n 8 -o out.vcf.gz input.vcf.gz \
  --- bcftools view -i "GT='alt'" -Ou \
  ::: bcftools view -s Sample -Oz
```

Chosen over `--` to avoid collision with tools (e.g. bcftools plugins) that use `--` as their own argument separator.

---

## Output file naming

### `-o` without `{}`

`{}` is absent — shard number is prepended as `shard_XX.`:

```
-o output.vcf.gz, -n 8   →  shard_01.output.vcf.gz ... shard_08.output.vcf.gz
-o out.bcf, -n 4          →  shard_1.out.bcf ... shard_4.out.bcf
-o results.txt, -n 100    →  shard_001.results.txt ... shard_100.results.txt
```

### `-o` with `{}`

`{}` is replaced with the zero-padded shard number:

```
-o output.{}.vcf.gz       →  output.01.vcf.gz ... output.08.vcf.gz
-o /results/{}/batch.bcf  →  /results/01/batch.bcf ... /results/08/batch.bcf
```

`{}` may appear anywhere in the path, including in a directory component. paravar creates all necessary parent directories (`mkdir -p`) before writing each shard.

### With `--gather`

`-o` is the final gather output path — no shard numbering applied. Temp shard files (written to `--tmp-dir`) are named using the `{}` or `shard_XX.` scheme derived from the `-o` value, with `.tmp` appended:

```
-o output.vcf.gz   →  $TMPDIR/paravar/paravar_shard_01.output.vcf.gz.tmp ...
-o out.{}.bcf      →  $TMPDIR/paravar/paravar_out.01.bcf.tmp ...
```

Temp files are deleted on success. On failure they are left on disk and their paths are printed to stderr.

### Format warning

If the input is VCF (`.vcf.gz`) and `-o` ends with `.bcf` (or vice versa), print a warning to stderr but proceed. Format conversion is the user's responsibility via the pipeline.

---

## Flags

### Common to all subcommands

| Flag | Long form | Description |
|---|---|---|
| `-o` | `--output` | Output path or suffix. Required for `scatter` and `run`. Optional for `run --gather` and `gather` — defaults to stdout if omitted. `/dev/stdout` also accepted. |
| `-v` | `--verbose` | Print progress to stderr |
| `-h` | `--help` | Show usage |

### `scatter` and `run`

| Flag | Long form | Description |
|---|---|---|
| `-n` | `--n-shards` | Number of shards (required, ≥ 1) |
| `-t` | `--max-threads` | Max threads for scatter/validation (default: min(n, 8)) |
| | `--force-scan` | Ignore index, scan all BGZF blocks — VCF only; exits 1 for BCF |

### `run` only

| Flag | Long form | Description |
|---|---|---|
| `-j` | `--max-jobs` | Max concurrent shard pipelines (default: n-shards) |
| | `--no-kill` | On failure, let sibling shards finish (default: kill siblings) |
| | `--gather` | Bare flag. Gather all shard outputs into single `-o` file |
| | `--tmp-dir <dir>` | Temp dir for gather (default: `$TMPDIR/paravar` or `/tmp/paravar`) |

### `run --gather` and `gather`

| Flag | Long form | Description |
|---|---|---|
| | `--header-pattern <pat>` | Strip lines with this prefix from shards 2..N — **text format only**; error if specified with VCF or BCF |
| | `--header-n <n>` | Strip first N lines from shards 2..N — **text format only**; error if specified with VCF or BCF |

`--header-pattern` and `--header-n` are mutually exclusive — exit 1 if both specified.

---

## Composability

All subcommands share the same output naming convention, so they compose naturally:

```bash
# Step by step
paravar scatter -n 8 -o out.vcf.gz input.vcf.gz
# → shard_01.out.vcf.gz ... shard_08.out.vcf.gz

paravar run -n 8 -o processed.vcf.gz input.vcf.gz --- bcftools view -i "GT='alt'" -Oz
# → shard_01.processed.vcf.gz ... shard_08.processed.vcf.gz

paravar gather -o merged.vcf.gz shard_*.processed.vcf.gz
# → merged.vcf.gz

# All in one
paravar run --gather -n 8 -o merged.vcf.gz input.vcf.gz --- bcftools view -i "GT='alt'" -Oz
# → merged.vcf.gz
```

---

## Gather logic (shared between `run --gather` and `gather`)

All gather behaviour lives in `gather.nim` and is shared between both modes.

### Format inference

From the `-o` extension (case-insensitive):

| Extension | Format | Compression |
|---|---|---|
| `.vcf.gz` or `.vcf.bgz` | VCF | BGZF |
| `.vcf` | VCF | Uncompressed |
| `.bcf` | BCF | BGZF |
| `.gz` or `.bgz` (any other stem) | Text | BGZF |
| Anything else | Text | Uncompressed |

Any extension not recognised as VCF or BCF is treated as text. `.gz` and `.bgz` are both treated as BGZF compression. No error is raised for unknown extensions.

### Format detection (sniffing)

Detected from the first bytes of the first shard's stream/file:

| First bytes (uncompressed) | Format |
|---|---|
| `BCF\x02\x02` | BCF |
| `##fileformat` | VCF |
| Anything else | Text |

For BGZF streams: decompress first block to check. Detected format stored globally; subsequent shards wait until detection completes.

### Header stripping

**VCF** (shards 2..N, automatic): strip all leading `#` lines.

**BCF** (shards 2..N, automatic): skip first `5 + 4 + l_text` uncompressed bytes (magic + l_text + header text). `--header-pattern`/`--header-n` are errors with VCF or BCF.

**Text** (shards 2..N, opt-in): `--header-pattern` or `--header-n`; default is no stripping.

### Header validation (`#CHROM` line check)

For VCF and BCF, before writing any output, paravar checks that the `#CHROM` line (the tab-separated sample column header) is byte-for-byte identical across all shards. This catches pipelines that reorder or add/remove samples.

**For `run --gather`**: the `#CHROM` line is extracted from shard 0's intercepted stream and stored. Each subsequent shard's interceptor extracts and compares before writing any data. On mismatch: exit 1 with a message identifying which shard failed and showing the differing lines. Kill all in-flight shard pipelines (respects `--no-kill`). Temp files left on disk.

**For `gather`**: extract the `#CHROM` line from each input file's header before concatenating anything. On first mismatch: exit 1 immediately with the differing shard path and lines. No partial output is written.

**For text format**: no header check — text gather is format-agnostic.

**VCF**: `#CHROM` line is the last `#`-prefixed line before data records. Extract from the decompressed header region.

**BCF**: `#CHROM` line is embedded in the `l_text` header blob. Extract by scanning the text for the line beginning with `#CHROM`.

### Recompression

| Incoming | Output compression | Action |
|---|---|---|
| Uncompressed | BGZF | Recompress via `compressToBgzfMulti` on the fly |
| BGZF | BGZF | Raw-copy blocks unchanged |
| BGZF | Uncompressed | Decompress, write raw bytes |
| Uncompressed | Uncompressed | Pass through unchanged |

### BGZF EOF handling

Each shard stream ends with a BGZF EOF block (28 bytes). Every interceptor/reader strips this trailing EOF before writing, so `concatenateShards` can write a single EOF at the end of the gather output. Temp shard files do not contain EOF blocks.

### Stdout output

When `-o` is omitted or is `/dev/stdout`, the gather output is written to stdout. Format inference still applies — infer from `/dev/stdout` is not meaningful, so when writing to stdout:
- Format and compression default to **uncompressed text** unless the sniffed shard format is VCF or BCF, in which case format matches the sniffed format with **no compression** (uncompressed VCF or uncompressed BCF stream)
- If the user wants a specific compressed format on stdout they should pipe: `paravar gather ... | bgzip -c > out.vcf.gz`

For `run --gather` writing to stdout: shard 0 writes directly to stdout fd. Shards 1..N still use temp files; `concatenateShards` writes them to stdout after all complete.

### Concatenation

1. Shard 0 → written directly to gather output or stdout (no temp file)
2. Shards 1..N → written to temp files
3. After all complete: open gather output in append mode (or write to stdout), raw-copy temp files in order
4. Write single BGZF EOF block if BGZF output (not applicable for stdout uncompressed)
5. Delete temp files on success; leave on disk and print paths to stderr on failure

For `gather` subcommand (operating on existing files): same logic, no temp files needed — read directly from input files.

---

## Module layout

| File | Responsibility |
|---|---|
| `src/paravar.nim` | Entry point |
| `src/paravar/main.nim` | CLI parsing, subcommand dispatch |
| `src/paravar/scatter.nim` | Scatter algorithm (VCF + BCF) |
| `src/paravar/bgzf_utils.nim` | Low-level BGZF I/O, no external deps beyond `-lz` |
| `src/paravar/run.nim` | Run mode: pipeline spawning, worker pool, interceptor coordination |
| `src/paravar/gather.nim` | Types, format inference, sniffing, stripping, interceptor thread, concatenation |

---

## Step-by-step plan

Execute in order. Do not skip ahead. Update checkboxes as steps complete.

### CLI refactor + gather (complete)
- [x] R0: reverted C1 partial implementation
- [x] C1–C6: subcommands, naming scheme, `:::` separator, `gather` subcommand, test updates

### Stdout + header validation (current milestone)
- [ ] **Step S1**: Implement stdout output for `run --gather` and `gather`: handle omitted `-o` and `/dev/stdout` as stdout fd. Format on stdout is uncompressed matching sniffed format. Write tests for stdout path (omitted `-o`, explicit `/dev/stdout`). Run `nimble test`.
- [ ] **Step S2**: Implement `#CHROM` header validation in `gather.nim`: `extractChromLine` for VCF and BCF, byte-for-byte comparison, early-exit on mismatch for both `run --gather` (interceptor) and `gather` (pre-scan all files before writing). Write tests covering match, mismatch (hard error), and text format (no check). Run `nimble test` — full suite green.

---

## Workflow rules for Claude Code

### Before starting any task

1. Re-read this file in full
2. Read the relevant source files
3. Consult `example/scatter_vcf.py` for scatter behaviour questions

### Hard rules

| Rule | Detail |
|---|---|
| **No new dependencies** | Do not add to `paravar.nimble` without asking |
| **Test before done** | Run `nimble test` and show full output before declaring any step complete |
| **No commits** | Stage changes, propose commit message, wait for user |
| **No layout changes** | Do not restructure modules without asking |
| **Small units** | One proc at a time, test it, proceed |
| **Ask when uncertain** | Ambiguous behaviour not covered here → stop and ask |

---

## Build reference

```bash
nimble build
nimble build -d:release
nimble test
nim c -d:debug -r tests/test_gather.nim
```

### BGZF EOF block constant

```nim
const BGZF_EOF* = [
  0x1f'u8, 0x8b, 0x08, 0x04, 0x00, 0x00, 0x00, 0x00,
  0x00, 0xff, 0x06, 0x00, 0x42, 0x43, 0x02, 0x00,
  0x1b, 0x00, 0x03, 0x00, 0x00, 0x00, 0x00, 0x00,
  0x00, 0x00, 0x00, 0x00
]
```

### BCF magic constant

```nim
const BCF_MAGIC* = [byte('B'), byte('C'), byte('F'), 0x02'u8, 0x02'u8]
```

---

## Out of scope (do not implement)

- Windows support
- bcftools as a gather dependency
- Tools that do not support stdin/stdout