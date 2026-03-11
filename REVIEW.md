# slide_snake Codebase Review

**Date**: 2026-03-11
**Scope**: Full pipeline — Snakefile, rules, Python/bash/R/awk scripts, configuration, documentation

---

## Summary of Completed Work

Four review passes fixed **~45 bugs**, implemented multiprocessing for barcode calling, reduced STAR I/O via `temp()` wrappers, enabled samtools threading across dedup rules, and corrected documentation mismatches. All fixes were committed in `e1ec3e6` plus pending uncommitted doc/script changes.

---

## Remaining Bugs

### MEDIUM-HIGH

| ID | File | Issue |
|----|------|-------|
| RB4 | `scripts/bash/fa_extract_rRNA.sh:77` | Unsafe `eval` on user-constructed command string — injection risk with special characters in paths. Replace with direct execution. |

### MEDIUM

| ID | File | Issue |
|----|------|-------|
| RB5 | `scripts/py/fastq_call_bc_umi_from_adapter_v2.py:226` | Hardcoded alignment score `420` when using integer adapter positions. Misleading for debugging. |
| RB6 | `scripts/bash/bam2splitBigWig.sh` | Missing `set -euo pipefail`. Piped commands (bedtools \| awk \| awk) fail silently. |
| RB7 | `scripts/bash/bam_chr2tag.sh:69-85` | No `set -eo pipefail` on `samtools view -h \| awk \| samtools view -bS` pipeline. |
| RB8 | `rules/short_read/5a_mirge.smk:61` | Only stderr captured. Add `1> {log.log}` for stdout. |
| RB9 | `scripts/bash/split_dedup.sh:135,142` | `$TEMP_LIST` unquoted in `samtools merge`. Use `"${TEMP_LIST[@]}"`. |

### LOW

| ID | File | Issue |
|----|------|-------|
| RB10 | `scripts/bash/bam_dedup.sh:83,90` | `[ -z ${INBAM} ]` should be `[ -z "${INBAM}" ]`. |
| RB11 | `rules/short_read/3a_star_align.smk:152-153` | Memory config type inconsistency between STAR passes. |
| RB12 | `scripts/py/fastq_internal_adapter_trim_R1_v2.py:149-151` | Entire read replaced with single `"N"` when adapter not found. |
| RB13 | Snakefile target list comprehensions | `RECIPE_DICT[SAMPLE]` without `.get()` fallback. |
| RB14 | `rules/0_utils.smk` | 15+ broad `except Exception:` blocks mask real errors. |
| RB15 | `scripts/py/cache_mtx_to_h5ad.py` | No error handling around file I/O; malformed input gives unhelpful errors. |

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
