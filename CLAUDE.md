# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

**slide_snake** is a Snakemake pipeline for spatial RNA-seq preprocessing and quantification. It supports multiple platforms (Visium, SlideSeq, StereoSeq, DBIT, DecoderSeq, miST) with both short-read (Illumina) and long-read (Oxford Nanopore) workflows.

## Common Commands

```bash
# Activate environment
mamba activate slsn

# Run pipeline locally
snakemake -k -p --use-conda --conda-frontend mamba -j 56

# Run on SLURM cluster
snakemake -k -p --use-conda --conda-frontend mamba --executor slurm \
  --workflow-profile profiles/slurm -j 24

# Run tests
python tests/test_barcode_processing.py

# Build documentation
sphinx-build -b html docs docs/_build/html
```

## Architecture

### Pipeline Flow

The Snakefile defines target outputs in `rule all`. Targets are organized by analysis type (barcode calling, alignment, QC) and are toggled by uncommenting lines. The pipeline has two parallel paths:

- **Short-read** (`rules/short_read/`): merge FASTQs → trim → call barcodes → rRNA filter → STAR/kallisto align → dedup → QC
- **Long-read** (`rules/ont/`): preprocess → trim → call barcodes → minimap2/ULTRA/isoquant align → QC

### Recipe System

The recipe system (`resources/recipe_sheet.csv`) is central to the pipeline's flexibility. Each recipe defines barcode/UMI extraction parameters, adapter sequences, and aligner-specific settings (STAR soloType, kallisto bus format, etc.). Recipes are looked up at runtime via `get_recipe_info(w, column, mode)` in `rules/0_utils.smk`. A single sample can use multiple recipes (comma-separated in the sample sheet).

### Key Files

- **Snakefile**: Orchestration, target definitions, rule includes
- **rules/0_utils.smk** (1200+ lines): All helper functions — sample sheet parsing, recipe lookups (`get_recipe_info`, `get_whitelist`, `get_fqs`, `get_bc_adapter`), validation (`check_sample_sheet`, `check_recipe_sheet`)
- **rules/0a_barcode_maps.smk**: Barcode whitelist processing rules
- **config.yaml**: Paths (sample sheet, output dir, tmp dir, recipe sheet) and global resource limits
- **profiles/default/config.yaml**: Per-rule thread/memory/runtime allocations (120+ rules specified)

### Barcode Processing Pipeline

1. **Extract** (`scripts/py/fastq_call_bc_umi_from_adapter_v2.py`): Uses parasail alignment to locate adapter sequences in reads and extract barcode + UMI
2. **Filter**: Remove reads where barcode extraction failed (marked with "-")
3. **Correct** (`scripts/py/tsv_bc_correction_parallelized.py`): Levenshtein distance-based correction against a whitelist, parallelized by kmer prefix
4. **Whitelist prep** (`scripts/py/process_barcode_whitelist.py`): Split/format barcodes for downstream tools (handles multi-component barcodes like SlideSeq 8bp+6bp)

### Sample Sheet

CSV with columns: `sampleID`, `fastq_R1`, `fastq_R2`, `recipe`, `ONT`, `recipe_ONT`, `BC_map`, `species`, `STAR_ref`, `genes_gtf`, `kb_idx`, `kb_t2g`, `genome_fa`, `cdna_fa`. See `docs/example_sample_sheet.csv`.

### Conda Environments

Each tool has its own env in `envs/`. The base environment is `envs/slsn.yml`. Snakemake's `--use-conda` handles activation per rule.

### Scripts Organization

- `scripts/py/`: Core processing (barcode calling, correction, matrix conversion, QC plotting)
- `scripts/bash/`: BAM manipulation, rRNA extraction, alignment wrappers
- `scripts/R/`: QC summary aggregation
- `scripts/awk/`: BAM/FASTQ field manipulation (tag extraction, filtering, trimming)

## Conventions

- Rule files are numbered by pipeline stage (1=preprocessing, 2=filtering/alignment, 3=post-alignment, 4+=quantification)
- Short-read and ONT rules are in separate subdirectories but share Python scripts via mode flags (e.g., mismatch thresholds differ)
- Wildcard `{SAMPLE}` is the primary sample identifier; `{RECIPE}` differentiates processing strategies per sample
- Output paths follow: `{OUTDIR}/{SAMPLE}/{tool}/{RECIPE}/...`
