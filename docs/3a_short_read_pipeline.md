# Short-Read Pipeline for `slide_snake`

The short-read pipeline in `slide_snake` is designed to process Illumina sequencing data from various spatial RNA-seq platforms. This pipeline handles the complete workflow from raw FASTQ files to count matrices and quality control reports.

## Overview

The short-read pipeline consists of several key stages:

1. **Preprocessing**: File merging, quality trimming, and barcode extraction
2. **rRNA Filtering**: Optional removal of ribosomal RNA contamination
3. **Alignment**: Mapping reads to reference genome using STAR or pseudoalignment with kallisto
4. **Quantification**: Generating gene expression count matrices
5. **Quality Control**: Comprehensive QC reporting at each step
6. **Advanced Analysis**: Optional small RNA analysis and uTAR detection

## Pipeline Architecture

```
Raw FASTQ files
    ↓
Preprocessing (merge, trim, extract barcodes)
    ↓
rRNA Filtering (optional)
    ↓
Alignment (STAR/kallisto)
    ↓
Quantification (gene counts)
    ↓
Quality Control & Visualization
    ↓
Count matrices & Reports
```

## Detailed Pipeline Steps

### Step 1: Preprocessing

#### 1a. Merge FASTQ Files
- **Purpose**: Combine multiple sequencing runs into single files
- **Input**: Multiple FASTQ files per read (R1, R2)
- **Output**: `merged_R1.fq.gz`, `merged_R2.fq.gz`
- **Features**:
  - Handles glob patterns for large numbers of files
  - Preserves read order and quality scores
  - Automatic compression detection and handling

#### 1b. Adapter Trimming
- **Tool**: `cutadapt`
- **Purpose**: Remove sequencing adapters and low-quality bases
- **Input**: Merged FASTQ files
- **Output**: `cut_R1.fq.gz`, `cut_R2.fq.gz`

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
- **Purpose**: Extract spatial barcodes and UMIs from reads
- **Method**: Recipe-specific strategies (fixed position, adapter-based, etc.)
- **Output**: `barcodes.tsv`, `barcodes_filtered.tsv`, `barcodes_corrected.tsv`

### Step 2: Ribosomal RNA Filtering (Optional)

#### 2a. BWA-based Filtering
- **Tool**: `BWA-MEM`
- **Process**:
  1. Align reads to rRNA reference
  2. Generate list of non-rRNA read IDs
  3. Filter R1/R2 files to remove rRNA reads
  4. Compress filtered files
- **Output**: `no_rRNA_R1.fq.gz`, `no_rRNA_R2.fq.gz`

#### 2b. RiboDetector Filtering
- **Tool**: `RiboDetector` (machine learning-based)
- **Advantages**: 
  - No reference required
  - Detects novel rRNA sequences
  - Higher sensitivity
- **Process**:
  1. Run RiboDetector on preprocessed reads
  2. Separate rRNA and non-rRNA reads
  3. Filter and compress output files

#### 2c. rRNA Alignment QC
- **Tool**: `Qualimap`
- **Purpose**: Assess rRNA contamination levels
- **Output**: `qualimap_rRNA_reports/`

### Step 3: STAR Alignment

#### 3a. STARsolo Alignment
- **Tool**: `STAR` with `--soloType` mode

- **Process**:
  1. Error correction of barcodes
  2. UMI counting
  3. Matrix generation
- **Output**: Count matrices in MTX format

### Step 5: Small RNA Analysis (Optional)

#### 5a. miRge3.0 Analysis
- **Tool**: `miRge3.0`
- **Purpose**: Quantify microRNAs and other small RNAs
- **Databases**:
  - **miRNAs**: miRBase v22 (hairpin & mature sequences)
  - **piRNAs**: piRBase v3.0 (gold standard sequences)
  - **Other**: tRNAs, snoRNAs, etc.

#### 5b. Custom Small RNA Pipeline
- **Tools**: `bowtie2` + custom scripts
- **Features**:
  - Hierarchical alignment to small RNA databases
  - Quantification with spatial barcode information
  - Integration with count matrices

### Step 6: uTAR Analysis (Unannotated Transcriptionally Active Regions)

#### 6a. uTAR Detection
- **Tool**: `uTAR_lite` (ported from [original](https://github.com/mckellardw/uTAR_lite))
- **Purpose**: Identify novel transcriptionally active regions
- **Method**: Uses `groHMM` to identify peaks in alignment data
- **Reference**: [Nature Communications 2021](https://www.nature.com/articles/s41467-021-22496-3)

#### 6b. uTAR Outputs
- **Annotations**: RefFlat and GTF files for detected uTARs
- **BAM files**: Tagged with uTAR assignments
- **Count matrices**: Both long-format TSV and sparse MTX formats

### Step 7: Data Integration and Caching

#### 7a. Seurat Integration
- **Tool**: `Seurat` (R)
- **Features**:
  - Spatial coordinates integration
  - Quality control metrics
  - Ready for downstream analysis
- **Output**: `.rds` files with spatial Seurat objects

#### 7b. Scanpy Integration
- **Tool**: `scanpy` (Python)
- **Features**:
  - AnnData object creation
  - Spatial coordinate handling
  - Quality control integration
- **Output**: `.h5ad` files ready for Python analysis

### Step 8: Quality Control and Reporting

#### 8a. Read-level QC
- **Tool**: `FastQC`
- **Stages**: 
  - Raw input (`preTrim`)
  - Post-trimming (`postTrim`)
  - Post-filtering (`postFilter`)
- **Custom**: Additional QC scripts for spatial-specific metrics

#### 8b. Alignment QC
- **Tool**: `Qualimap`
- **Features**:
  - RNA-seq specific metrics
  - Spatial barcode statistics
  - Coverage analysis

#### 8c. Unmapped Read Analysis
- **FastQC**: Quality assessment of unmapped reads
- **BLAST**: Alignment of most abundant unmapped sequences
- **Purpose**: Identify contamination or missing references

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
out/
├── {SAMPLE}/
│   ├── illumina/
│   │   ├── fastq/
│   │   │   ├── merged/
│   │   │   └── cut/
│   │   ├── STAR/
│   │   │   └── {RECIPE}/
│   │   │       ├── Solo.out/
│   │   │       └── Aligned.sortedByCoord.out.bam
│   │   ├── kb/
│   │   │   └── {RECIPE}/
│   │   │       └── counts_unfiltered/
│   │   └── qc/
│   │       ├── fastqc/
│   │       └── qualimap/
│   └── bc/
│       └── whitelist.txt
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
- BED and GTF files for TARs.
- Tagged and sorted BAM files.
- Count matrices in both long-format (`counts.tsv.gz`) and sparse matrix format (`uTAR.mtx.gz`).


### Step 4: Pseudoalignment w/ kallisto/bustools
#TODO


### Step 5: Small/Micro RNA Analysis
- Align reads to small RNA databases using `miRge3.0` and generate count matrices.
- Output: `miRge3_reports`, `count_matrices`

**Work-in-progress...**
- Currently have rules set up to align to databases listed below (with `bowtie2`) and generate a count matrix (note, this uses `STAR`'s cell/bead/spot barcoding)
- Need to better optimize alignment params

#### miRge3.0
- [Link to documentation](https://mirge3.readthedocs.io/en/latest/quick_start.html)
- [Link to library download](https://sourceforge.net/projects/mirge3/files/miRge3_Lib/)
- Bulk miRg3.0 is implemented, but it is quite slow. I recommend commenting out the target rule(s) unless you have small/total RNA data for which you need this analysis.

#### Small RNA Reference Databases
- **piRNAs** - [pirbase](http://bigdata.ibp.ac.cn/piRBase/), v3.0
  - *Note* - use "gold standard" piRs because there are way too many predicted sequences...
- **miRNAs** - [mirbase](https://mirbase.org/), v22
  - Use both hairpin & mature sequences, in that order


### Step 6: Seurat & scanpy caching
- Use `Seurat` and `scanpy` to generate a cache of the count matrices for downstream analysis.
- Objects have spatial coordinates and are ready for downstream analysis.


### Step 7: Quality Control 
- Perform quality control on the reads and alignments using `fastqc` and custom scripts (`readqc`).
- Output: `fastqc_reports`

