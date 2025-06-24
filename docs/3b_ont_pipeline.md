
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
ONT FASTQ files
    ↓
Preprocessing (merge, filter, strand)
    ↓
Barcode and UMI Calling
    ↓
Read Trimming and Quality Control
    ↓
Alignment (Genome or Transcriptome)
    ↓
Quantification and Visualization
    ↓
Count matrices and Reports
```

## Detailed Pipeline Steps

### Step 1: Preprocessing (`rules/ont/1a_preprocessing.smk`)

#### 1a. File Format Merging
- **Rule**: `ont_1a_merge_formats`
- **Purpose**: Consolidate multiple ONT output files into single FASTQ
- **Features**:
  - Handles different ONT output formats
  - Preserves quality scores and read headers
  - Automatic compression detection

#### 1b. Length Filtering  
- **Rule**: `ont_1a_length_filter`
- **Purpose**: Remove reads below minimum length threshold
- **Parameters**: Configurable minimum read length (typically 200-500bp)
- **Rationale**: Short reads often lack sufficient information for barcode calling

#### 1c. Adapter Scanning and Stranding
- **Rule**: `ont_1a_call_adapter_scan`  
- **Purpose**: Identify adapter sequences and determine read orientation
- **Method**: Searches for forward and reverse primer sequences
- **Output**: Stranded FASTQ files (`forward.fq.gz`, `reverse.fq.gz`)

**Key Adapter Sequences:**
- **Forward primer**: Near barcode region (e.g., `CTACACGACGCTCTTCCGATCT`)
- **Reverse primer**: Template switch oligo (e.g., `AAGCAGTGGTATCAACGCAGAG`)

#### 1d. Read Splitting by Poly(T) Stretches
- **Rule**: `ont_1a_split_fastq_to_R1_R2`
- **Purpose**: Separate reads into barcode (R1) and transcript (R2) components
- **Method**: Identifies poly(T) sequences that separate barcode and cDNA regions
- **Output**: `R1.fq.gz` (barcodes/UMIs), `R2.fq.gz` (transcript sequences)

### Step 2: Trimming (`rules/ont/1b_trimming.smk`)

#### 2a. Cutadapt Trimming
- **Rule**: `ont_1b_cutadapt`
- **Purpose**: Remove adapter sequences and low-quality regions
- **Features**:
  - Platform-specific adapter sequences
  - Quality score trimming
  - Minimum length filtering

#### 2b. Hard Trimming (R1)
- **Rule**: `ont_1b_R1_hardTrimming`
- **Purpose**: Remove fixed-length sequences from barcode reads
- **Use case**: When adapter positions are consistent

#### 2c. Internal Adapter Trimming
- **Rule**: `ont_1b_R1_internalTrim`
- **Purpose**: Remove internal adapter sequences
- **Method**: Searches for and removes adapter sequences within reads

#### 2d. Trimming Summary
- **Rule**: `ont_1b_cutadapt_summary`
- **Output**: Visualization of trimming statistics
- **Metrics**: Adapter removal rates, length distributions, quality improvements

### Step 3: Barcode and UMI Calling (`rules/ont/1c_barcode_calling.smk`)

This is the most complex and critical step for ONT spatial data processing.

#### 3a. Adapter-Based Barcode Calling
- **Rule**: `ont_1c_fastq_call_bc_from_adapter`
- **Method**: Locates adapter sequences and extracts barcodes/UMIs relative to adapter positions
- **Features**:
  - Multiple barcode components (e.g., combinatorial indexing)
  - Flexible positioning (left or right of adapter)
  - Configurable edit distance tolerance

**Configuration Parameters:**
```yaml
BC_adapter: "TCTTCAGCGTTCCCGAGA TCTTCAGCGTTCCCGAGA"  # Adapter sequences
BC_start: "0 26"                                      # Start positions
BC_length: "8 6"                                      # Barcode lengths  
BC_offset: "0 0"                                      # Offsets from adapter
BC_position: "left right"                            # Position relative to adapter
BC_max_ED: 2                                         # Max edit distance
BC_min_ED_diff: 1                                    # Min difference between top matches
```

#### 3b. Barcode Filtering
- **Rule**: `ont_1c_filter_barcodes`
- **Purpose**: Remove low-quality barcode calls
- **Criteria**:
  - Minimum quality scores
  - Edit distance thresholds
  - Adapter match confidence

#### 3c. Barcode Correction
- **Rule**: `ont_1c_tsv_bc_correction`
- **Method**: Correct barcodes against whitelist using edit distance
- **Features**:
  - Configurable correction stringency
  - Multi-component barcode handling
  - Concatenation strategies for complex barcodes

#### 3d. Correction Summary
- **Rule**: `ont_1c_summarize_bc_correction`
- **Output**: Statistics on barcode correction success rates
- **Metrics**: 
  - Correction rates per barcode component
  - Edit distance distributions
  - Whitelist match statistics

### Step 4: Alignment

The ONT pipeline supports both genome and transcriptome alignment strategies.

#### 4a. Genome Alignment (`rules/ont/1d_minimap2_genome.smk`)

##### Junction Preparation
- **Rule**: `ont_1d_genome_generate_junction_bed`
- **Purpose**: Convert GTF annotations to junction BED format
- **Use**: Improves splice-aware alignment accuracy

##### Minimap2 Genome Alignment
- **Rule**: `ont_1d_genome_align_minimap2_genome`
- **Tool**: `minimap2` with `-ax splice` mode
- **Features**:
  - Splice-aware alignment
  - Long-read optimized parameters
  - Customizable alignment settings

##### Post-processing
- **Rule**: `ont_1d_genome_sort_compress_output`
- **Purpose**: Sort and compress alignment outputs
- **Rule**: `ont_1d_genome_add_corrected_barcodes`
- **Purpose**: Add corrected barcode information to BAM tags

#### 4b. Transcriptome Alignment (`rules/ont/1d_minimap2_transcriptome.smk`)

##### Minimap2 Transcriptome Alignment
- **Rule**: `ont_1d_txome_align_minimap2_transcriptome`
- **Tool**: `minimap2` with `-ax map-ont` mode
- **Purpose**: Direct alignment to transcript sequences
- **Advantages**: 
  - Faster than genome alignment
  - Direct transcript assignment
  - Suitable for transcript quantification

##### Barcode Sorting
- **Rule**: `ont_1d_txome_sort_by_cb`
- **Purpose**: Sort alignments by cell barcode for efficient processing

##### Oarfish Quantification
- **Rule**: `ont_1d_txome_oarfish_quant`
- **Tool**: `oarfish` (from COMBINE-lab)
- **Purpose**: Transcript-level quantification from long reads
- **Features**:
  - EM algorithm for multi-mapping reads
  - UMI-aware quantification
  - Transcript isoform resolution

#### 4c. Kallisto Long-Read Processing (`rules/ont/1e_kallisto-lr.smk`)

##### Kallisto-lr Quantification
- **Rule**: `ont_1e_kallisto_lr`
- **Tool**: `kallisto` with long-read mode
- **Purpose**: Fast pseudoalignment and quantification
- **Output**: BUS format files for downstream processing

##### BUS to Matrix Conversion
- **Rule**: `ont_1e_bus2mat_lr`
- **Tool**: `bustools`
- **Purpose**: Convert BUS files to count matrices
- **Features**:
  - Gene-level aggregation
  - UMI counting
  - Sparse matrix output

### Step 5: Advanced Quantification (`rules/ont/1g_isoquant.smk`)

#### IsoQuant Analysis
- **Rule**: `ont_1g_isoquant`
- **Tool**: `IsoQuant`
- **Purpose**: Comprehensive transcript and gene quantification
- **Features**:
  - Novel isoform detection
  - Gene expression quantification
  - Alternative splicing analysis

#### Gene Assignment
- **Rule**: `ont_1g_add_isoquant_genes_to_bam`
- **Purpose**: Add gene assignment tags to BAM files
- **Use**: Enables gene-level filtering and analysis

#### UMI-tools Counting
- **Rule**: `ont_1g_umitools_count`
- **Tool**: `UMI-tools`
- **Purpose**: Generate final count matrices
- **Features**:
  - UMI deduplication
  - Spatial barcode integration
  - Multiple output formats

### Step 6: Quality Control (`rules/ont/2a_readqc.smk`, `rules/ont/2b_qualimap.smk`)

#### Read-Level QC
- **Raw input QC**: `ont_2a_readQC_0_rawInput`
- **Pre-trimming QC**: `ont_2a_readQC_1_preCutadapt`  
- **Post-trimming QC**: `ont_2a_readQC_2_postCutadapt`

**QC Metrics:**
- Read length distributions
- Quality score profiles
- Adapter content
- Base composition

#### Alignment QC
- **Qualimap RNA-seq**: `ont_2b_qualimap` - RNA-seq specific metrics
- **Qualimap BAM QC**: `ont_2b_qualimap_bamqc` - General alignment statistics

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
BC_start: 0
BC_length: 16
UMI_length: 12
```

### Complex Barcode System (Multi-component)
```yaml
BC_adapter: "ADAPTER1 ADAPTER2"
BC_start: "0 26"
BC_length: "8 6"
BC_position: "left right"
BC_concat: True
```

## Output Structure

```
out/
├── {SAMPLE}/
│   ├── ont/
│   │   ├── preprocessing/
│   │   │   ├── merged.fq.gz
│   │   │   ├── length_filtered.fq.gz
│   │   │   └── stranded/
│   │   ├── barcodes/
│   │   │   ├── raw_barcodes.tsv
│   │   │   ├── filtered_barcodes.tsv
│   │   │   └── corrected_barcodes.tsv
│   │   ├── minimap2_genome/
│   │   │   └── {RECIPE}/
│   │   │       └── aligned.sorted.bam
│   │   ├── minimap2_transcriptome/
│   │   │   └── {RECIPE}/
│   │   │       ├── aligned.sorted.bam
│   │   │       └── oarfish_quant/
│   │   ├── kb_lr/
│   │   │   └── {RECIPE}/
│   │   │       └── counts/
│   │   └── isoquant/
│   │       └── {RECIPE}/
│   │           ├── counts.tsv
│   │           └── novel_isoforms.gtf
│   └── qc/
│       ├── ont_readqc/
│       └── qualimap/
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

This comprehensive ONT pipeline provides robust processing of long-read spatial RNA-seq data with flexible barcode handling and multiple quantification strategies.

