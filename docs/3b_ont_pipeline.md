
## Table of Contents
1. [Preprocessing](#preprocessing)
2. [Barcode and UMI Calling](#barcode-and-umi-calling)
3. [Trimming](#trimming)
4. [Alignment](#alignment)
   - [Genome Alignment](#genome-alignment)
   - [Transcriptome Alignment](#transcriptome-alignment)
5. [Quantification](#quantification)
6. [Quality Control](#quality-control)
7. [Visualization](#visualization)

---

## Preprocessing

### File: `rules/ont/1a_preprocessing.smk`

- **Purpose**: Prepares raw ONT reads for downstream analysis.
- **Key Rules**:
  - `ont_1a_merge_formats`: Merges input files into a single `.fastq.gz` file.
  - `ont_1a_length_filter`: Filters out reads shorter than a specified minimum length.
  - `ont_1a_call_adapter_scan`: Identifies adapter sequences and generates stranded `.fastq.gz` files.
  - `ont_1a_split_fastq_to_R1_R2`: Splits reads into R1 and R2 based on poly(T) stretches.


## Trimming

### File: `rules/ont/1b_trimming.smk`

- **Purpose**: Trims adapters and low-quality regions from ONT reads.
- **Key Rules**:
  - `ont_1b_cutadapt`: Trims adapters and low-quality bases using Cutadapt.
  - `ont_1b_R1_hardTrimming`: Performs hard trimming on R1 reads.
  - `ont_1b_R1_internalTrim`: Removes internal adapter sequences from R1 reads.
  - `ont_1b_cutadapt_summary`: Generates a summary plot of trimming results.


## Barcode and UMI Calling

### File: `rules/ont/1c_barcode_calling.smk`

- **Purpose**: Extracts and corrects barcodes and UMIs from ONT reads.
- **Key Rules**:
  - `ont_1c_fastq_call_bc_from_adapter`: Calls barcodes and UMIs based on adapter sequences.
  - `ont_1c_filter_barcodes`: Filters out invalid barcodes.
  - `ont_1c_tsv_bc_correction`: Corrects barcodes using a whitelist.
  - `ont_1c_summarize_bc_correction`: Summarizes barcode correction statistics.


## Alignment

### Genome Alignment

#### File: `rules/ont/1d_minimap2_genome.smk`

- **Purpose**: Aligns ONT reads to a genome reference using Minimap2.
- **Key Rules**:
  - `ont_1d_genome_generate_junction_bed`: Converts GTF to junction BED for spliced alignment.
  - `ont_1d_genome_align_minimap2_genome`: Aligns reads to the genome.
  - `ont_1d_genome_sort_compress_output`: Sorts and compresses alignment output.
  - `ont_1d_genome_add_corrected_barcodes`: Adds corrected barcodes to BAM files.

### Transcriptome Alignment

#### File: `rules/ont/1d_minimap2_transcriptome.smk`

- **Purpose**: Aligns ONT reads to a transcriptome reference using Minimap2.
- **Key Rules**:
  - `ont_1d_txome_align_minimap2_transcriptome`: Aligns reads to the transcriptome.
  - `ont_1d_txome_sort_by_cb`: Sorts BAM files by cell barcode.
  - `ont_1d_txome_oarfish_quant`: Quantifies transcripts using Oarfish.


#### File: `rules/ont/1e_kallisto-lr.smk`
- **Purpose**: Processes and visualizes long-read data using Kallisto.
- **Key Rules**:
  - `ont_1e_kallisto_lr`: Runs Kallisto for long-read quantification.


## Quantification

### File: `rules/ont/1g_isoquant.smk`

- **Purpose**: Performs transcript and gene quantification using IsoQuant.
- **Key Rules**:
  - `ont_1g_isoquant`: Runs IsoQuant for transcript and gene quantification.
  - `ont_1g_add_isoquant_genes_to_bam`: Adds gene tags to BAM files.
  - `ont_1g_umitools_count`: Generates count matrices using UMI-tools.


## Quality Control

### File: `rules/ont/2a_readqc.smk`

- **Purpose**: Performs quality control on raw and processed reads.
- **Key Rules**:
  - `ont_2a_readQC_0_rawInput`: QC on raw input reads.
  - `ont_2a_readQC_1_preCutadapt`: QC before adapter trimming.
  - `ont_2a_readQC_2_postCutadapt`: QC after adapter trimming.

### File: `rules/ont/2b_qualimap.smk`

- **Purpose**: Performs alignment QC using Qualimap.
- **Key Rules**:
  - `ont_2b_qualimap`: Runs Qualimap RNA-seq QC.
  - `ont_2b_qualimap_bamqc`: Runs Qualimap BAM QC.

