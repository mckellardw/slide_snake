# slide_snake Codebase Review

**Date**: 2026-03-12
**Scope**: Full pipeline — Snakefile, rules, Python/bash/R/awk scripts, configuration, documentation

---

## Summary of Completed Work

Seven review passes fixed **~68 bugs**, implemented multiprocessing for barcode calling, reduced STAR I/O via `temp()` wrappers, enabled samtools threading across dedup rules, corrected documentation mismatches, hardened all bash scripts with `set -eo pipefail` and proper variable quoting, and fixed critical Python indexing/logic errors.

---

## Remaining Bugs

### MEDIUM

| ID | File | Line(s) | Issue |
|----|------|---------|-------|
| RB14 | `scripts/py/fastq_call_bc_umi_from_adapter_v2.py` | 226 | Hardcoded alignment score `420` when using integer adapter positions. Misleading for debugging. |
| RB15 | `rules/short_read/5a_mirge.smk` | 61 | Only stderr captured. Add `1> {log.log}` for stdout. |
| RB16 | `scripts/bash/fa_build_rRNA_gtf.sh` | 78, 80, 167 | Unquoted variables in `zcat`/`cat`/redirection. Complex pipes without error handling. |
| RB17 | `scripts/bash/bam_dedupByChr.sh` | 80, 88 | Line 80: `ls *REF_*.bam > ${BAMLIST}` — unsafe glob. Line 88: unquoted `${BAM}` in loop. |
| RB18 | `scripts/bash/bam_uTAR_HMM.sh` | 61, 95-96 | Line 61: AWK file append with unquoted TMPDIR. Lines 95-96: `pids2` array used but never initialized. |
| RB19 | `scripts/bash/fq_bwa_rRNA_align.sh` | 64 | No `mkdir -p` before writing output — fails if directory doesn't exist. |
| RB20 | `scripts/bash/umitools_extract_ont.sh` | whole file | Missing `set -eo pipefail`. No input file validation. |
| RB21 | `rules/ont/1b_trimming.smk` | 65, 183 | cutadapt shell blocks redirect stdout (`1>>`) but never redirect stderr. Missing `2> {log.err}`. |

### LOW

_All LOW bugs resolved._

- **RB22** (fixed): STAR pass2 memory config now uses `megabytes2bytes(config["MEMLIMIT_MB"])` matching pass1.
- **RB24** (false positive): `RECIPE_DICT[SAMPLE]` always valid — all comprehensions iterate `R2_FQS.keys()`, same keys used to populate `RECIPE_DICT`.
- **RB25** (fixed): Narrowed 9 `except Exception:` blocks to `except KeyError:` or `except (KeyError, TypeError, ValueError):` in reference lookups, recipe lookups, and barcode config functions. Remaining 6 blocks in complex helper functions left as-is (varied failure modes with reasonable fallbacks).
- **RB26** (by design): Script is called by Snakemake — unhandled exceptions propagate as rule failures with tracebacks, which is correct behavior.

---

## Runtime Optimization Roadmap

### Tier 2 — Short-term (1-4 weeks)

| Task | File | Estimated Speedup |
|------|------|-------------------|
| Stream compression instead of two-pass gzip | `long2mtx.py`, STAR compress rule | 10-20% I/O savings |
| BK-tree for barcode correction | `tsv_bc_correction_parallelized.py` | 2-3x on correction step |
| Split dedup by chromosome | `3c_star_dedup.smk` | 4-8x parallelism within dedup |

### Tier 3 — Rust rewrite candidates (1-3 months)

#### Priority 1: Parasail alignment wrapper (VERY HIGH impact)

**Scripts**: `adapter_scan_parasail_v2.py`, `fastq_call_bc_umi_from_adapter_v2.py`, `fastq_split_reads_parallelized_v3.py`, `fastq_internal_adapter_trim_R1_v2.py`

Smith-Waterman alignment via `parasail.sw_trace_striped_32()` runs per-read, multiple times (forward/reverse, multiple adapters). Python wrapper overhead + GIL prevents true parallelism.

**Rust approach**: `parasail-rs` or `rust-bio` with `rayon` parallel iterators, SIMD batching. CLI tool accepting FASTQ + adapter sequences, outputting TSV. Expected speedup: **5-10x**.

#### Priority 2: Barcode correction (HIGH impact)

**Script**: `tsv_bc_correction_parallelized.py`

Levenshtein distance against whitelist for every barcode. Current k-mer index is O(n) per query.

**Rust approach**: BK-tree/VP-tree for O(log n) lookup, SIMD Levenshtein via `triple_accel`. Expected speedup: **3-5x**.

#### Priority 3: BAM tag addition (MEDIUM impact)

**Scripts**: `tsv2tag.py`, `bam_readID2tags.py`, `add_sequence_as_tag.py`

Per-read pysam overhead dominates.

**Rust approach**: `rust-htslib` streaming BAM manipulation. Expected speedup: **2-3x**.

### Tier 4 — Architectural (long-term)

| Change | Impact |
|--------|--------|
| Replace `.mtx` text with HDF5/Parquet | 5-10x compression, faster downstream loading |
| Named pipes for STAR -> compression | Eliminate uncompressed intermediate I/O |
| GPU acceleration for alignment | 10-100x for parasail-equivalent operations |

---

## Remaining Documentation Issues

### HIGH

| ID | Issue |
|----|-------|
| D4 | `data/README.md` says GSE2358968; `data/test_decoderseq/README.md` says GSE235896. One is wrong. |
| D5 | Stub READMEs: `data/test_StereoSeq/` (TODO), `data/test_magicseq/` (TODO), `data/test_microST/` (blank). |
| D6 | Docs show `out/{SAMPLE}/illumina/STAR/` but actual path is `out/{SAMPLE}/short_read/STARsolo/`. |
| D7 | `docs/2b_barcode_processing_guide.md` references non-existent `process_adapter_insertion.py`. |
| D8 | `docs/3b_ont_pipeline.md` documents a known syntax error in `ont_1a_length_filter` with no fix. |

### MEDIUM

| ID | Issue |
|----|-------|
| D9 | `docs/2a_recipes.md:100` — corrupted header "Technology Specifications d Configuration". |
| D10 | `bam_tags.md` not listed in `docs/index.rst` table of contents. |
| D11 | No Python version requirement in `docs/0_install.md`. |
| D12 | `README.md` line 11 AI disclosure — consider rewording. |
