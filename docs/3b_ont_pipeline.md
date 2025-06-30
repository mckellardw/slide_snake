# Oxford Nanopore Technologies (ONT) Pipeline for `slide_snake`

The ONT pipeline in `slide_snake` processes long-read sequencing data from Oxford Nanopore Technologies platforms for spatial RNA-seq applications. This pipeline is specifically designed to handle the unique challenges of ONT data including variable read lengths, higher error rates, and complex barcode structures.

## Overview

The ONT pipeline addresses several key challenges specific to long-read spatial RNA-seq:
- **Variable barcode positions** due to synthesis artifacts
- **Complex adapter structures** requiring flexible extraction strategies  
- **Higher error rates** necessitating robust correction algorithms
- **Full-length transcript detection** enabling isoform-level analysis
- **UMI integration** for accurate quantification

## Pipeline Architecture

```
ONT BAM/FASTQ files
    ↓
Step 1: Preprocessing (merge, stranding, read-splitting, & trimming)
    ↓
Step 2: Barcode and UMI Calling
    ↓
Step 3: Alignment (Genome or Transcriptome) & quantification
    ↓
Step 4: Scanpy/Seurat objects and QC
```

## Detailed Pipeline Steps

### Step 1: Preprocessing (`rules/ont/1a_preprocessing.smk`)

#### 1a. File Format Merging
- **Rule**: `ont_1a_merge_formats`
- **Purpose**: Consolidate multiple ONT output files into single FASTQ
- **Input**: Multiple ONT files (various formats)
- **Output**: `merged.fq.gz`
- **Features**:
  - Handles different ONT output formats
  - Preserves quality scores and read headers
  - Automatic compression detection

#### 1c. Adapter Scanning and Stranding
- **Rule**: `ont_1a_call_adapter_scan`  
- **Purpose**: Identify adapter sequences and determine read orientation
- **Input**: `merged.fq.gz`
- **Output**: 
  - `adapter_scan.tsv` (adapter detection results)
  - `merged_stranded.fq.gz` (stranded reads)
- **Method**: Searches for forward and reverse primer sequences
- **Tools**: Uses parasail for alignment-based adapter detection

**Key Adapter Sequences:**
- **Forward primer**: Near barcode region (e.g., `CTACACGACGCTCTTCCGATCT`)
- **Reverse primer**: Template switch oligo (e.g., `AAGCAGTGGTATCAACGCAGAG`)

#### 1d. Read Splitting by Poly(T) Stretches
- **Rule**: `ont_1a_split_fastq_to_R1_R2`
- **Purpose**: Separate reads into barcode (R1) and transcript (R2) components
- **Input**: `merged_adapter.fq.gz`
- **Output**: 
  - `merged_adapter_R1.fq.gz` (barcodes/UMIs)
  - `merged_adapter_R2.fq.gz` (transcript sequences)
  - `merged_adapter_ambiguous.fq.gz` (ambiguous splits)
- **Method**: Identifies poly(T) sequences that separate barcode and cDNA regions
- **Parameters**: Uses 8bp poly(T) sequence as default split anchor

#### Additional Processing Steps

The pipeline includes several additional intermediate steps:

- **Read ID filtering**: `ont_1a_readIDs_by_adapter_type` creates lists of read IDs for different adapter types
- **Read list merging**: `ont_1a_merge_scan_lists` combines read ID lists for downstream processing  
- **FASTQ subsetting**: `ont_1a_subset_fastq_by_adapter_type` filters reads based on adapter scan results
- **Compression**: `ont_1a_compress_merged_fq` compresses intermediate FASTQ files
- **BAM filtering**: Multiple steps filter BAM files to remove reads without proper barcodes and UMIs

### Step 2: Trimming (`rules/ont/1b_trimming.smk`)

#### 2a. Cutadapt Trimming
- **Rule**: `ont_1b_cutadapt`
- **Purpose**: Remove adapter sequences and low-quality regions
- **Input**: 
  - `merged_adapter_R1.fq.gz`
  - `merged_adapter_R2.fq.gz`
- **Output**: 
  - `cut_R1.fq.gz`
  - `cut_R2.fq.gz`
  - `1b_cutadapt.json` (trimming statistics)
- **Features**:
  - Platform-specific adapter sequences
  - Quality score trimming (minimum Q5 for ONT)
  - Minimum length filtering
  - TSO removal and polyadenylation handling

#### 2b. Hard Trimming (R1)
- **Rule**: `ont_1b_R1_hardTrimming`
- **Purpose**: Remove fixed-length sequences from barcode reads
- **Input**: `cut_R1.fq.gz`
- **Output**: `cut_hardTrim_R1.fq.gz`
- **Use case**: When adapter positions are consistent
- **Parameters**: Fixed positions for CB1 and CB2 boundaries

#### 2c. Internal Adapter Trimming
- **Rule**: `ont_1b_R1_internalTrim`
- **Purpose**: Remove internal adapter sequences
- **Input**: `cut_R1.fq.gz`
- **Output**: `cut_internalTrim_R1.fq.gz`
- **Method**: Searches for and removes adapter sequences within reads
- **Features**: Uses parasail for alignment-based trimming

#### 2d. Trimming Summary
- **Rule**: `ont_1b_cutadapt_summary`
- **Input**: `1b_cutadapt.json`
- **Output**: `1b_cutadapt_summary.png` (visualization of trimming statistics)
- **Metrics**: Adapter removal rates, length distributions, quality improvements

### Step 3: Barcode and UMI Calling (`rules/ont/1c_barcode_calling.smk`)

This is the most complex and critical step for ONT spatial data processing.

#### 3a. Adapter-Based Barcode Calling
- **Rule**: `ont_1c_fastq_call_bc_from_adapter`
- **Input**: Trimmed R1 FASTQ files
- **Output**: 
  - `barcodes.tsv` (raw barcode calls)
  - `barcode_stats.tsv` (calling statistics)
- **Method**: Locates adapter sequences and extracts barcodes/UMIs relative to adapter positions
- **Features**:
  - Multiple barcode components (e.g., combinatorial indexing)
  - Flexible positioning (left or right of adapter)
  - Configurable edit distance tolerance (default: 0.7 similarity)

**Configuration Parameters:**
```yaml
BC_adapter: "TCTTCAGCGTTCCCGAGA"              # Adapter sequence
BC_length: "8"                                # Barcode length  
BC_offset: "0"                                # Offset from adapter
BC_position: "left"                           # Position relative to adapter
BC_max_ED: 2                                  # Max edit distance
BC_min_ED_diff: 1                            # Min difference between top matches
```

#### 3b. Barcode Filtering
- **Rule**: `ont_1c_filter_barcodes`
- **Purpose**: Remove low-quality barcode calls
- **Input**: `barcodes.tsv`
- **Output**: `barcodes_filtered.tsv`
- **Criteria**: Removes entries with missing barcode information ("-" entries)

#### 3c. Barcode Correction
- **Rule**: `ont_1c_tsv_bc_correction`
- **Input**: 
  - `barcodes_filtered.tsv`
  - Barcode whitelist files
- **Output**: 
  - `barcodes_corrected.tsv` (corrected barcodes only)
  - `barcodes_corrected_full.tsv` (detailed correction info)
- **Method**: Correct barcodes against whitelist using edit distance
- **Features**:
  - Configurable correction stringency
  - Multi-component barcode handling
  - Concatenation strategies for complex barcodes

#### 3d. Correction Summary
- **Rule**: `ont_1c_summarize_bc_correction`
- **Input**: `barcodes_corrected_full.tsv`
- **Output**: `bc_correction_stats.txt` (statistics on barcode correction success rates)
- **Metrics**: 
  - Correction rates per barcode component
  - Edit distance distributions
  - Whitelist match statistics

### Step 4: Alignment and Quantification

The ONT pipeline supports both genome and transcriptome alignment strategies.

#### 4a. Genome Alignment (`rules/ont/2a_minimap2_genome.smk`)

##### Junction Preparation
- **Rule**: `ont_2a_genome_generate_junction_bed`
- **Purpose**: Convert GTF annotations to junction BED format
- **Input**: GTF annotation file
- **Output**: `junctions.bed`
- **Use**: Improves splice-aware alignment accuracy

##### Minimap2 Genome Alignment
- **Rule**: `ont_2a_genome_align_minimap2_genome`
- **Tool**: `minimap2` with `-ax splice` mode
- **Input**: 
  - R2 FASTQ file (transcript sequences)
  - `junctions.bed`
- **Output**: `tmp.sam` (temporary SAM file)
- **Features**:
  - Splice-aware alignment
  - Long-read optimized parameters
  - Customizable alignment settings

##### Post-processing
- **Rule**: `ont_2a_genome_sort_compress_output`
- **Purpose**: Sort and compress alignment outputs
- **Input**: `tmp.sam`
- **Output**: `sorted.bam`

- **Rule**: `ont_2a_genome_add_corrected_barcodes`
- **Purpose**: Add corrected barcode information to BAM tags
- **Input**: 
  - `sorted.bam`
  - `barcodes_corrected.tsv`
- **Output**: `sorted_cb.bam`

#### 4b. Transcriptome Alignment (`rules/ont/2b_minimap2_transcriptome.smk`)

##### Minimap2 Transcriptome Alignment
- **Rule**: `ont_2b_txome_align_minimap2_transcriptome`
- **Tool**: `minimap2` with `-ax map-ont` mode
- **Input**: R2 FASTQ file (transcript sequences)
- **Output**: `aligned.bam`
- **Purpose**: Direct alignment to transcript sequences
- **Advantages**: 
  - Faster than genome alignment
  - Direct transcript assignment
  - Suitable for transcript quantification

##### Barcode and UMI Addition
- **Rule**: `ont_2b_txome_add_corrected_barcodes`
- **Input**: 
  - `aligned.bam`
  - `barcodes_corrected.tsv`
- **Output**: `aligned_cb.bam`

- **Rule**: `ont_2b_txome_add_umis`
- **Input**: 
  - `aligned_cb.bam`
  - `barcodes_filtered.tsv`
- **Output**: `aligned_cb_ub.bam`

##### Barcode Sorting
- **Rule**: `ont_2b_txome_sort_by_cb`
- **Purpose**: Sort alignments by cell barcode for efficient processing

##### Oarfish Quantification
- **Rule**: `ont_2b_txome_oarfish_quant`
- **Tool**: `oarfish` (from COMBINE-lab)
- **Purpose**: Transcript-level quantification from long reads
- **Features**:
  - EM algorithm for multi-mapping reads
  - UMI-aware quantification
  - Transcript isoform resolution

#### 4c. Kallisto Long-Read Processing (`rules/ont/2c_kallisto-lr.smk`)

##### Kallisto-lr Quantification
- **Rule**: `ont_2c_kallisto_lr`
- **Tool**: `kallisto` with long-read mode
- **Purpose**: Fast pseudoalignment and quantification
- **Output**: BUS format files for downstream processing

##### BUS to Matrix Conversion
- **Rule**: `ont_2c_bus2mat_lr`
- **Tool**: `bustools`
- **Purpose**: Convert BUS files to count matrices
- **Features**:
  - Gene-level aggregation
  - UMI counting
  - Sparse matrix output

### Step 5: Advanced Quantification (`rules/ont/2e_isoquant.smk`)

#### IsoQuant Analysis
- **Rule**: `ont_2e_isoquant`
- **Tool**: `IsoQuant`
- **Purpose**: Comprehensive transcript and gene quantification
- **Features**:
  - Novel isoform detection
  - Gene expression quantification
  - Alternative splicing analysis

#### UMI-tools Counting
- **Rule**: `ont_2e_umitools_count`
- **Tool**: `UMI-tools`
- **Input**: Final processed BAM file with all tags
- **Output**: `umitools_counts.tsv.gz` (UMI-deduplicated count matrix)
- **Purpose**: Generate final count matrices
- **Features**:
  - UMI deduplication
  - Spatial barcode integration
  - Gene-level quantification

#### Matrix Format Conversion
- **Rule**: `ont_2e_counts_to_sparse`
- **Input**: `umitools_counts.tsv.gz`
- **Output**: 
  - `barcodes.tsv.gz`
  - `features.tsv.gz`  
  - `matrix.mtx.gz`
- **Purpose**: Convert long-format counts to sparse matrix format

#### Gene Assignment
- **Rule**: `ont_2e_add_isoquant_genes_to_bam`
- **Purpose**: Add gene assignment tags to BAM files
- **Input**: 
  - Filtered BAM file
  - IsoQuant read assignments
- **Output**: BAM with gene tags (IG)
- **Use**: Enables gene-level filtering and analysis

### Step 6: Quality Control (`rules/ont/3a_readqc.smk`, `rules/ont/3b_qualimap.smk`)

#### Read-Level QC
- **Raw input QC**: `ont_3a_readQC_0_rawInput`
  - **Input**: `merged.fq.gz`
  - **Output**: `0_rawInput/merged_qc.tsv`
- **Pre-trimming QC**: `ont_3a_readQC_1_preCutadapt`
  - **Input**: `merged_adapter_R1.fq.gz`, `merged_adapter_R2.fq.gz`
  - **Output**: `1_preCutadapt/{READ}_qc.tsv`
- **Post-trimming QC**: `ont_3a_readQC_2_postCutadapt`
  - **Input**: `cut_R1.fq.gz`, `cut_R2.fq.gz`
  - **Output**: `2_postCutadapt/{READ}_qc.tsv`

**QC Metrics:**
- Read length distributions
- Quality score profiles
- Adapter content
- Base composition

#### Alignment QC
- **Qualimap RNA-seq**: `ont_3b_qualimap` - RNA-seq specific metrics
- **Qualimap BAM QC**: `ont_3b_qualimap_bamqc` - General alignment statistics

**Alignment Metrics:**
- Mapping rates
- Insert size distributions  
- Gene coverage profiles
- Transcript assignment rates

## ONT-Specific Considerations

### Error Rate Management
- **Barcode calling**: Uses edit distance tolerance and confidence metrics
- **Alignment**: Optimized parameters for higher error rates
- **Quantification**: UMI-based correction for technical noise

### Adapter Complexity
- **Multiple adapters**: Handles complex adapter structures with multiple components
- **Variable positions**: Adapter-based extraction accommodates synthesis variations
- **Quality filtering**: Removes low-confidence adapter matches

### Full-Length Transcript Analysis
- **Isoform detection**: Capable of detecting novel splice variants
- **Alternative splicing**: Quantifies different transcript isoforms
- **Long-range connectivity**: Links distant exons within single reads

## Configuration Examples

### Standard Visium ONT
```yaml
recipe_ONT: "visium_ont"
fwd_primer: "CTACACGACGCTCTTCCGATCT"
rev_primer: "AAGCAGTGGTATCAACGCAGAG"
BC_adapter: "TCTTCAGCGTTCCCGAGA"
BC_length: 16
BC_offset: 0
BC_position: "left"
BC_max_ED: 2
BC_min_ED_diff: 1
UMI_adapter: "TCTTCAGCGTTCCCGAGA"
UMI_length: 12
UMI_offset: 16
UMI_position: "left"
```

### Complex Barcode System (Multi-component)
```yaml
BC_adapter: "ADAPTER1 ADAPTER2"
BC_length: "8 6"
BC_offset: "0 0"
BC_position: "left right"
BC_max_ED: "2 2"
BC_min_ED_diff: "1 1"
BC_concat: True
```

## Output Structure

```
out/
├── {SAMPLE}/
│   ├── ont/
│   │   ├── tmp/                              # Temporary processing files
│   │   │   ├── merged.fq.gz
│   │   │   ├── merged_adapter_R1.fq.gz
│   │   │   ├── merged_adapter_R2.fq.gz
│   │   │   ├── cut_R1.fq.gz
│   │   │   └── cut_R2.fq.gz
│   │   ├── adapter_scan.tsv                  # Adapter detection results
│   │   ├── barcodes_umis/
│   │   │   └── {RECIPE}/
│   │   │       ├── barcodes.tsv              # Raw barcode calls
│   │   │       ├── barcodes_filtered.tsv     # Filtered barcodes
│   │   │       ├── barcodes_corrected.tsv    # Corrected barcodes
│   │   │       └── bc_correction_stats.txt   # Correction statistics
│   │   ├── minimap2/
│   │   │   └── {RECIPE}/
│   │   │       ├── junctions.bed
│   │   │       └── sorted_cb.bam             # Genome-aligned BAM with barcodes
│   │   ├── minimap2_txome/
│   │   │   └── {RECIPE}/
│   │   │       ├── aligned_cb_ub.bam         # Transcriptome-aligned BAM
│   │   │       └── oarfish_quant/            # Oarfish quantification results
│   │   ├── kb_lr/
│   │   │   └── {RECIPE}/
│   │   │       └── counts/                   # Kallisto long-read counts
│   │   ├── isoquant/
│   │   │   └── {RECIPE}/
│   │   │       ├── counts.tsv
│   │   │       └── novel_isoforms.gtf
│   │   ├── readqc/                           # Quality control results
│   │   │   ├── 0_rawInput/
│   │   │   ├── 1_preCutadapt/
│   │   │   └── 2_postCutadapt/
│   │   ├── plots/                            # Visualization outputs
│   │   │   ├── 1a_adapter_scan_summary.png
│   │   │   └── 1b_cutadapt_summary.png
│   │   └── logs/                             # Processing logs
│   └── qc/
│       └── qualimap/                         # Alignment quality metrics
```

## Best Practices

### Recipe Selection
- **Start simple**: Begin with standard recipes before trying complex configurations
- **Platform matching**: Use recipes designed for your specific platform
- **Quality assessment**: Always check QC reports before downstream analysis

### Performance Optimization
- **Memory management**: ONT processing can be memory-intensive
- **Parallel processing**: Most steps support multi-threading
- **Storage**: Plan for large intermediate files

### Quality Control
- **Barcode correction rates**: Aim for >60% correction success
- **Mapping rates**: Should be >70% for good quality data
- **UMI complexity**: Check for appropriate UMI diversity

### Troubleshooting
- **Low barcode detection**: Adjust adapter sequences and edit distance parameters
- **Poor mapping**: Verify reference transcriptome completeness
- **Memory errors**: Reduce parallel jobs or increase memory allocation

## Known Issues and Notes

### Current Implementation Notes
- The `ont_1a_length_filter` rule has a syntax error with duplicate input directives
- Some intermediate steps create temporary files that are automatically cleaned up
- The pipeline supports multiple alignment strategies (genome vs transcriptome) that can be run in parallel
- Memory requirements can be substantial, particularly for large datasets during alignment and barcode correction

### Pipeline Flexibility
- Most rules support recipe-specific configurations through the recipe sheet
- Barcode calling parameters can be tuned per-recipe for different spatial platforms
- Multiple quantification strategies are available (IsoQuant, Oarfish, Kallisto-lr) depending on analysis goals

This comprehensive ONT pipeline provides robust processing of long-read spatial RNA-seq data with flexible barcode handling and multiple quantification strategies.

