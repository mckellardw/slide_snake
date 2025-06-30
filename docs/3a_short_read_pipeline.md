# Short-Read Pipeline for `slide_snake`

The short-read pipeline in `slide_snake` is designed to process Illumina sequencing data from various spatial RNA-seq platforms. This pipeline handles the complete workflow from raw FASTQ files to count matrices and quality control reports.

## Overview

The short-read pipeline consists of several key stages:

1. **Preprocessing** (`1a-1c`): File merging, quality trimming, and barcode extraction
2. **rRNA Filtering** (`2a-2c`, optional): Removal of ribosomal RNA contamination
3. **Alignment & Quantification** (`3a-3u`): Mapping reads to reference genome using STAR or pseudoalignment with kallisto
4. **Pseudoalignment** (`4a`): Alternative quantification using kallisto/bustools
5. **Small RNA Analysis** (`5a`, optional): microRNA and small RNA quantification
6. **Final Outputs** (`6a-6b`): Converting to standard formats (Scanpy/Seurat)
7. **Quality Control** (`7a-7b`): Comprehensive QC reporting at each step
8. **Advanced Analysis**: Optional uTAR detection and unmapped read analysis

## Pipeline Architecture

```
Raw FASTQ files
    ↓
Preprocessing (merge, trim, extract barcodes)
    ↓
rRNA Filtering (optional)
    ↓
Alignment & quantification (STAR/kallisto)
    ↓
Scanpy/Seurat objects and QC
```

## Detailed Pipeline Steps

### Step 1: Preprocessing

#### 1a. Merge FASTQ Files
- **Rule**: `ilmn_1a_merge_fastqs`
- **Purpose**: Combine multiple sequencing runs into single files
- **Input**: Multiple FASTQ files per read (R1, R2) from sample sheet
- **Output**: 
  - `{OUTDIR}/{SAMPLE}/short_read/tmp/merged_R1.fq.gz` (merged R1 reads)
  - `{OUTDIR}/{SAMPLE}/short_read/tmp/merged_R2.fq.gz` (merged R2 reads)
- **Features**:
  - Handles glob patterns for large numbers of files
  - Preserves read order and quality scores
  - Automatic compression detection and handling
  - Uses temporary files to save disk space

#### 1b. Adapter Trimming
- **Rule**: `ilmn_1b_cutadapt`
- **Tool**: `cutadapt`
- **Purpose**: Remove sequencing adapters and low-quality bases
- **Input**: 
  - `{OUTDIR}/{SAMPLE}/short_read/tmp/merged_R1.fq.gz`
  - `{OUTDIR}/{SAMPLE}/short_read/tmp/merged_R2.fq.gz`
- **Output**: 
  - `{OUTDIR}/{SAMPLE}/short_read/tmp/cut_R1.fq.gz` (trimmed R1 reads)
  - `{OUTDIR}/{SAMPLE}/short_read/tmp/cut_R2.fq.gz` (trimmed R2 reads)
  - `{OUTDIR}/{SAMPLE}/short_read/logs/cutadapt1.json` (trimming statistics)

**Adapter Sequences Removed:**
- **5' Adapters**:
  - SlideSeq TSO reverse complement: `CCCATTCACTCTGCGTTGATACCAGCTT`
- **3' Adapters**:
  - Poly-A homopolymers: `100×A`
  - Poly-G homopolymers: `100×G` (important for 2-color sequencers)
  - Poly-T homopolymers: `100×T`
  - Nextera sequences: `CTGTCTCTTATA`, `TATAAGAGACAG`
  - Curio TSO: `AAGCTGGTATCAACGCAGAGTGAATGGG`
  - Curio R1 adapter: `TCTTCAGCGTTCCCGAGA`
  - Illumina universal: `AGATCGGAAGAG`

#### 1c. Custom Barcode and UMI Extraction
- **Rule**: `ilmn_1c_fastq_call_bc_from_adapter`
- **Purpose**: Extract spatial barcodes and UMIs from reads
- **Input**: Trimmed FASTQ files
- **Output**: 
  - `{OUTDIR}/{SAMPLE}/bc/barcodes.tsv` (raw barcode calls)
  - `{OUTDIR}/{SAMPLE}/bc/barcode_stats.tsv` (extraction statistics)
- **Method**: Recipe-specific strategies (fixed position, adapter-based, etc.)

**Additional Barcode Processing:**
- **Rule**: `ilmn_1c_filter_barcodes`
  - **Input**: `{OUTDIR}/{SAMPLE}/bc/barcodes.tsv`
  - **Output**: `{OUTDIR}/{SAMPLE}/bc/barcodes_filtered.tsv` (removes missing barcode entries)

- **Rule**: `ilmn_1c_tsv_bc_correction`
  - **Input**: `{OUTDIR}/{SAMPLE}/bc/barcodes_filtered.tsv` + whitelist files
  - **Output**: 
    - `{OUTDIR}/{SAMPLE}/bc/barcodes_corrected.tsv` (corrected barcodes)
    - `{OUTDIR}/{SAMPLE}/bc/barcodes_corrected_full.tsv` (detailed correction info)

### Step 2: Ribosomal RNA Filtering (Optional)

#### 2a. BWA-based Filtering
- **Tool**: `BWA-MEM`
- **Rules**: 
  - `ilmn_2a_extract_rRNA_fasta` - Extract rRNA sequences from reference
  - `ilmn_2a_build_rRNA_bwa_index` - Build BWA index for rRNA
  - `ilmn_2a_align_rRNA_bwa` - Align reads to rRNA reference
  - `ilmn_2a_get_rRNA_unaligned_reads` - Extract non-rRNA reads
- **Process**:
  1. Extract rRNA sequences from cDNA FASTA using keywords
  2. Build BWA index for rRNA sequences
  3. Align reads to rRNA reference
  4. Generate list of non-rRNA read IDs
  5. Filter R1/R2 files to remove rRNA reads
- **Input**: `{OUTDIR}/{SAMPLE}/short_read/tmp/cut_R1.fq.gz`, `{OUTDIR}/{SAMPLE}/short_read/tmp/cut_R2.fq.gz`
- **Output**: 
  - `{OUTDIR}/{SAMPLE}/short_read/tmp/no_rRNA_R1.fq.gz` (rRNA-depleted R1 reads)
  - `{OUTDIR}/{SAMPLE}/short_read/tmp/no_rRNA_R2.fq.gz` (rRNA-depleted R2 reads)
  - `{OUTDIR}/{SAMPLE}/short_read/bwa_rRNA/rRNA_sequences.fa.gz` (extracted rRNA reference)

#### 2b. RiboDetector Filtering  
- **Rule**: `ilmn_2b_ribodetector`
- **Tool**: `RiboDetector` (machine learning-based)
- **Advantages**: 
  - No reference required
  - Detects novel rRNA sequences
  - Higher sensitivity
- **Input**: `{OUTDIR}/{SAMPLE}/short_read/tmp/cut_R1.fq.gz`, `{OUTDIR}/{SAMPLE}/short_read/tmp/cut_R2.fq.gz`
- **Output**: 
  - `{OUTDIR}/{SAMPLE}/short_read/tmp/cut_R1_nonrRNA.fq.gz` (filtered R1)
  - `{OUTDIR}/{SAMPLE}/short_read/tmp/cut_R2_nonrRNA.fq.gz` (filtered R2)
- **Process**:
  1. Run RiboDetector on preprocessed reads
  2. Separate rRNA and non-rRNA reads
  3. Filter and compress output files

#### 2c. rRNA Alignment QC
- **Rule**: `ilmn_2c_rRNA_qualimap`
- **Tool**: `Qualimap`
- **Purpose**: Assess rRNA contamination levels
- **Input**: rRNA-aligned BAM files
- **Output**: `{OUTDIR}/{SAMPLE}/short_read/qc/qualimap_rRNA/` (QC metrics for rRNA content)

### Step 3: STAR Alignment

#### 3a. STARsolo Alignment (Standard Pipeline)
- **Rule**: `ilmn_3a_STARsolo_align`
- **Tool**: `STAR` with `--soloType` mode
- **Input**: Processed FASTQ files + barcode whitelists
- **Output**:
  - `{OUTDIR}/{SAMPLE}/short_read/STARsolo/{RECIPE}/Aligned.sortedByCoord.out.bam` (coordinate-sorted alignments)
  - `{OUTDIR}/{SAMPLE}/short_read/STARsolo/{RECIPE}/Solo.out/Gene/raw/matrix.mtx` (gene count matrix)
  - `{OUTDIR}/{SAMPLE}/short_read/STARsolo/{RECIPE}/Solo.out/Gene/raw/barcodes.tsv` (cell barcodes)
  - `{OUTDIR}/{SAMPLE}/short_read/STARsolo/{RECIPE}/Solo.out/Gene/raw/features.tsv` (gene features)
  - `{OUTDIR}/{SAMPLE}/short_read/STARsolo/{RECIPE}/Solo.out/GeneFull/raw/` (full gene counts including introns)
  - `{OUTDIR}/{SAMPLE}/short_read/STARsolo/{RECIPE}/Solo.out/Velocyto/raw/` (spliced/unspliced matrices for RNA velocity)
  - `{OUTDIR}/{SAMPLE}/short_read/STARsolo/{RECIPE}/Unmapped.out.mate1`, `Unmapped.out.mate2` (unmapped reads)
- **Features**:
  - Cell barcode error correction
  - UMI counting and deduplication
  - Multiple count matrix types (Gene, GeneFull, Velocyto)

#### 3b. Custom Barcode Pipeline
- **Rule**: `ilmn_3e_STARcustom_firstPass` (2-pass alignment, first pass)
- **Rule**: `ilmn_3e_STARcustom_secondPass` (2-pass alignment, second pass)
- **Purpose**: For custom barcode extraction workflows
- **Input**: Processed FASTQ files
- **Output**: 
  - `{OUTDIR}/{SAMPLE}/short_read/STAR_custom/{RECIPE}/Aligned.sortedByCoord.out.bam` (alignments without barcode tags)
  - `{OUTDIR}/{SAMPLE}/short_read/STAR_custom/{RECIPE}/firstPass/_STARpass1/SJ.out.tab` (splice junctions from first pass)

#### 3c. Barcode Integration
- **Rule**: `ilmn_3e_STARcustom_add_barcode_tags`
- **Purpose**: Add custom barcode calls to BAM files
- **Input**: 
  - `{OUTDIR}/{SAMPLE}/short_read/STAR_custom/{RECIPE}/Aligned.sortedByCoord.out.bam`
  - `{OUTDIR}/{SAMPLE}/bc/barcodes_corrected.tsv`
- **Output**: `{OUTDIR}/{SAMPLE}/short_read/STAR_custom/{RECIPE}/Aligned.sortedByCoord.tagged.bam` (BAM with barcode tags)

#### 3d. UMI-tools Quantification
- **Rule**: `ilmn_3e_STARcustom_umitools_count`
- **Tool**: `UMI-tools`
- **Purpose**: Generate count matrices from tagged BAM files
- **Input**: Tagged BAM file
- **Output**: `{OUTDIR}/{SAMPLE}/short_read/STAR_custom/{RECIPE}/umitools_counts.tsv.gz` (UMI-deduplicated counts)

### Step 5: Small RNA Analysis (Optional)

#### 5a. miRge3.0 Analysis
- **Rule**: `ilmn_5a_mirge3`
- **Tool**: `miRge3.0`
- **Purpose**: Quantify microRNAs and other small RNAs
- **Input**: Processed FASTQ files
- **Output**:
  - `{OUTDIR}/{SAMPLE}/short_read/miRge3/{RECIPE}/miRge.*.csv` (miRNA quantification)
  - `{OUTDIR}/{SAMPLE}/short_read/miRge3/{RECIPE}/annotation.report.html` (detailed report)
- **Databases**:
  - **miRNAs**: miRBase v22 (hairpin & mature sequences)
  - **piRNAs**: piRBase v3.0 (gold standard sequences)
  - **Other**: tRNAs, snoRNAs, etc.
- **Features**:
  - Hierarchical alignment to small RNA databases
  - Isotype classification and novel miRNA prediction
  - Comprehensive annotation reports

### Step 6: uTAR Analysis (Unannotated Transcriptionally Active Regions)

#### 6a. uTAR Detection
- **Rule**: `ilmn_3u_STAR_uTAR` 
- **Tool**: `uTAR_lite` (ported from [original](https://github.com/mckellardw/uTAR_lite))
- **Purpose**: Identify novel transcriptionally active regions
- **Method**: Uses `groHMM` to identify peaks in alignment data
- **Input**: STAR-aligned BAM files
- **Output**:
  - `{OUTDIR}/{SAMPLE}/short_read/uTAR/{RECIPE}/uTAR_refFlat.txt` (RefFlat format annotations)
  - `{OUTDIR}/{SAMPLE}/short_read/uTAR/{RECIPE}/uTAR.gtf` (GTF format annotations)
  - `{OUTDIR}/{SAMPLE}/short_read/uTAR/{RECIPE}/Aligned.sortedByCoord.tagged.bam` (BAM with uTAR tags)
  - `{OUTDIR}/{SAMPLE}/short_read/uTAR/{RECIPE}/counts.tsv.gz` (long-format counts)
  - `{OUTDIR}/{SAMPLE}/short_read/uTAR/{RECIPE}/uTAR.mtx.gz` (sparse matrix format)
- **Reference**: [Nature Communications 2021](https://www.nature.com/articles/s41467-021-22496-3)

### Step 7: Data Integration and Caching

#### 7a. Scanpy Integration
- **Rule**: `ilmn_6a_cache_h5ad_STAR` (for STAR outputs)
- **Rule**: `ilmn_6a_cache_h5ad_KB` (for kallisto outputs)
- **Tool**: `scanpy` (Python)
- **Input**: Count matrices (MTX format) + barcode/feature files
- **Output**: 
  - `{OUTDIR}/{SAMPLE}/short_read/STARsolo/{RECIPE}/Solo.out/{SOLO}/raw/{ALGO}.h5ad` (AnnData objects)
  - `{OUTDIR}/{SAMPLE}/short_read/STARsolo/{RECIPE}/Solo.out/{SOLO}/raw/{ALGO}_qc_plots.png` (QC plots)
- **Features**:
  - AnnData object creation with spatial coordinates
  - Quality control metrics integration
  - Ready for Python-based downstream analysis

#### 7b. Seurat Integration
- **Rule**: `ilmn_6b_cache_rds_STAR` (for STAR outputs)
- **Rule**: `ilmn_6b_cache_rds_KB` (for kallisto outputs)
- **Tool**: `Seurat` (R)
- **Input**: Count matrices + spatial coordinate files
- **Output**: 
  - `{OUTDIR}/{SAMPLE}/short_read/STARsolo/{RECIPE}/Solo.out/{SOLO}/raw/{ALGO}.rds` (Seurat objects)
  - `{OUTDIR}/{SAMPLE}/short_read/STARsolo/{RECIPE}/Solo.out/{SOLO}/raw/{ALGO}_qc_plots.png` (QC plots)
- **Features**:
  - Spatial Seurat object creation
  - Integrated spatial coordinates and metadata
  - Quality control metrics
  - Ready for R-based downstream analysis

### Step 8: Quality Control and Reporting

#### 8a. Read-level QC
- **Rules**: 
  - `ilmn_7a_fastqc_preTrim` (raw input QC)
  - `ilmn_7a_fastqc_postTrim` (post-trimming QC)
  - `ilmn_7a_fastqc_postFilter` (post-filtering QC)
  - `ilmn_7a_fastqc_unmapped` (unmapped reads QC)
- **Tool**: `FastQC`
- **Input**: FASTQ files at various processing stages
- **Output**: 
  - `{OUTDIR}/{SAMPLE}/short_read/qc/fastqc_preTrim_reports/` (raw data QC)
  - `{OUTDIR}/{SAMPLE}/short_read/qc/fastqc_postTrim_reports/` (post-trimming QC)
  - `{OUTDIR}/{SAMPLE}/short_read/qc/fastqc_postFilter_reports/` (post-filtering QC)
  - `{OUTDIR}/{SAMPLE}/short_read/qc/fastqc_unmapped_reports/` (unmapped reads QC)
- **Custom QC**: 
  - **Rule**: `ilmn_7b_readqc`
  - **Output**: `{OUTDIR}/{SAMPLE}/short_read/qc/readqc_reports/` (spatial-specific metrics)

#### 8b. Alignment QC
- **Rule**: `ilmn_3q_qualimap_genome`
- **Tool**: `Qualimap`
- **Input**: Aligned BAM files
- **Output**: `{OUTDIR}/{SAMPLE}/short_read/qc/qualimap_genome/` (alignment QC reports)
- **Features**:
  - RNA-seq specific metrics
  - Spatial barcode statistics
  - Coverage analysis
  - Insert size distributions

#### 8c. Unmapped Read Analysis
- **Rules**: 
  - `ilmn_7a_fastqc_unmapped` (FastQC on unmapped reads)
  - `ilmn_3b_STAR_unmapped_BLAST` (BLAST analysis of unmapped reads)
- **Purpose**: Identify contamination or missing references
- **Output**:
  - `{OUTDIR}/{SAMPLE}/short_read/qc/fastqc_unmapped_reports/` (unmapped read quality)
  - `{OUTDIR}/{SAMPLE}/short_read/STAR_unmapped/{RECIPE}/blast_results/` (contamination analysis)

## Configuration and Recipes

### Recipe Selection
Different recipes optimize for different platforms and data types:

```yaml
# Standard mRNA analysis
recipe: "visium"

# Total RNA analysis (includes non-coding RNAs)
recipe: "visium_total"

# rRNA filtering
recipe: "visium_std_rRNA-bwa"

# Benchmarking multiple approaches
recipe: "visium visium_total seeker_MatchLinker"
```

### Memory and Performance Optimization

#### Resource Configuration
```yaml
# config.yaml
CORES: 16              # CPU cores per job
MEMLIMIT_GB: 32G       # Memory limit for STAR
MEMLIMIT_MB: 32000     # Memory limit (MB) for other tools
```

#### Platform-Specific Optimizations
- **STAR**: Optimized for splice-aware alignment
- **Kallisto**: Fast pseudoalignment, good for quantification
- **Memory management**: Automatic temporary file cleanup

## Output Directory Structure

```
{OUTDIR}/
├── {SAMPLE}/
│   ├── short_read/
│   │   ├── tmp/                              # Temporary files (auto-cleaned)
│   │   │   ├── merged_R1.fq.gz              # Merged R1 reads
│   │   │   ├── merged_R2.fq.gz              # Merged R2 reads
│   │   │   ├── cut_R1.fq.gz                 # Trimmed R1 reads
│   │   │   ├── cut_R2.fq.gz                 # Trimmed R2 reads
│   │   │   ├── no_rRNA_R1.fq.gz             # rRNA-filtered R1 (if enabled)
│   │   │   └── no_rRNA_R2.fq.gz             # rRNA-filtered R2 (if enabled)
│   │   ├── STARsolo/
│   │   │   └── {RECIPE}/
│   │   │       ├── Aligned.sortedByCoord.out.bam
│   │   │       ├── Solo.out/
│   │   │       │   ├── Gene/raw/            # Gene-level counts
│   │   │       │   ├── GeneFull/raw/        # Gene + intron counts
│   │   │       │   └── Velocyto/raw/        # Spliced/unspliced counts
│   │   │       ├── Unmapped.out.mate1
│   │   │       └── Unmapped.out.mate2
│   │   ├── STAR_custom/
│   │   │   └── {RECIPE}/
│   │   │       ├── Aligned.sortedByCoord.out.bam
│   │   │       ├── Aligned.sortedByCoord.tagged.bam
│   │   │       └── umitools_counts.tsv.gz
│   │   ├── kbpython_std/
│   │   │   └── {RECIPE}/
│   │   │       └── counts_unfiltered/
│   │   │           └── cellranger/         # kallisto outputs
│   │   ├── miRge3/
│   │   │   └── {RECIPE}/
│   │   │       ├── miRge.*.csv             # miRNA quantification
│   │   │       └── annotation.report.html
│   │   ├── uTAR/
│   │   │   └── {RECIPE}/
│   │   │       ├── uTAR_refFlat.txt
│   │   │       ├── uTAR.gtf
│   │   │       ├── Aligned.sortedByCoord.tagged.bam
│   │   │       ├── counts.tsv.gz
│   │   │       └── uTAR.mtx.gz
│   │   ├── qc/
│   │   │   ├── fastqc_preTrim_reports/
│   │   │   ├── fastqc_postTrim_reports/
│   │   │   ├── fastqc_postFilter_reports/
│   │   │   ├── fastqc_unmapped_reports/
│   │   │   ├── qualimap_genome/
│   │   │   ├── qualimap_rRNA/
│   │   │   └── readqc_reports/
│   │   ├── bwa_rRNA/
│   │   │   └── rRNA_sequences.fa.gz        # Extracted rRNA reference
│   │   └── logs/                           # Log files for all steps
│   └── bc/
│       ├── barcodes.tsv                    # Raw barcode calls
│       ├── barcodes_filtered.tsv           # Filtered barcodes
│       ├── barcodes_corrected.tsv          # Error-corrected barcodes
│       ├── barcodes_corrected_full.tsv     # Detailed correction info
│       ├── barcode_stats.tsv               # Extraction statistics
│       └── whitelist.txt                   # Platform-specific whitelist
```

## Best Practices

### For New Users
1. **Start with test data**: Use provided examples to validate installation
2. **Choose appropriate recipe**: Match recipe to your platform
3. **Monitor resource usage**: Adjust memory/CPU based on your system
4. **Validate outputs**: Check QC reports before downstream analysis

### For Production Analysis
1. **Use rRNA filtering**: Especially for total RNA libraries
2. **Enable multiple recipes**: For benchmarking and validation
3. **Set up proper references**: Use high-quality, version-controlled references
4. **Archive parameters**: Keep record of exact recipes and versions used

### Troubleshooting
- **Low cell detection**: Check barcode whitelist and recipe parameters
- **High rRNA contamination**: Enable rRNA filtering steps
- **Memory errors**: Reduce parallel jobs or increase memory limits
- **Adapter issues**: Verify adapter sequences in recipe configuration

This comprehensive short-read pipeline provides robust, reproducible processing of spatial RNA-seq data with extensive quality control and flexibility for different platforms and experimental designs.

