# Short read preprocessing info



### kallisto/bustools
- Barcode & UMI parameters for `kallisto bus`:
```bash
-x 0,0,14:0,14,21:1,0,0
```


## QC
- `fastqc` is run on R2 files before (`preTrim`) and after (`postTrim`) adapter trimming
- `qualimap` is used to assess `STAR` alignment

## Unmapped Read Analysis
- `fastqc` is run on all unmapped reads from `STAR`
- The most abundant unmapped reads are also aligned with `blast`


## Pipeline Steps

### Step 1: Preprocessing
#### (a) Merge Fastq Files
- Merge multiple sequencing runs into a single .fastq file for each read.
- Can specify multiple input files for each read, using regular expressions for large numbers of files.
- Output: `merged_R1.fq.gz`, `merged_R2.fq.gz`

#### (b) Trimming
- Trim adapters and low-quality bases from the reads using `cutadapt`.
- Output: `cut_R1.fq.gz`, `cut_R2.fq.gz`

**Five Prime Adapter(s):**
- Reverse complement of the SlideSeq TSO to reduce strand invasion artifacts: `CCCATTCACTCTGCGTTGATACCAGCTT`

**Three Prime Adapter(s):**
- "A" homopolymers: `100-A`
- "G" homopolymers, important for Illumina 2-color sequencers: `100-G`
- "T" homopolymers: `100-T`
- Nextera adapter sequence: `CTGTCTCTTATA`
- Reverse complement of Nextera sequence: `TATAAGAGACAG`
- Curio template switch oligo (TSO) - remove any polyadenylated TSOs: `AAGCTGGTATCAACGCAGAGTGAATGGG`
- Curio R1 internal adapter - shows up in R2, and used for R1 trimming: `TCTTCAGCGTTCCCGAGA`
- Reverse of Curio R1 adapter: `AGAGCCCTTGCGACTTCT`
- Illumina universal sequence: `AGATCGGAAGAG`

#### (c) Custom Barcode and UMI Extraction
**work-in-progress...**
- Extract barcodes and UMIs from the reads.
- Output: `barcodes.tsv`, `barcodes_filtered.tsv`, `barcodes_corrected.tsv`

### Step 2: Ribosomal RNA Filtering
- Filter out ribosomal RNA reads using either `RiboDetector` or `BWA`.
- Output: `no_rRNA_R1.fq.gz`, `no_rRNA_R2.fq.gz`

#### (a) BWA
- BWA is used to align reads to the rRNA reference and filter out rRNA reads.
- The following steps are performed:
  - Align reads to the rRNA reference using BWA.
  - Generate a list of read IDs to keep (non-rRNA) from the filtered R2 file.
  - Filter R1 reads using the generated list of read IDs.
  - Compress the filtered fastq files.
  - 
#### (b) RiboDetector
- RiboDetector is used to filter out ribosomal RNA reads from the preprocessed reads.
- The following steps are performed:
  - Run RiboDetector on preprocessed reads to separate rRNA and non-rRNA reads.
  - Generate a list of read IDs to keep (non-rRNA) from the output R2 file.
  - Filter R1 reads using the generated list of read IDs.
  - Compress the filtered fastq files.


#### (c) Qualimap for rRNA Alignment
- Use `qualimap` to assess the alignment of rRNA reads to the reference genome.
- Output: `qualimap_rRNA_reports`

### Step 3: STAR Alignment
- Align reads to the reference genome using `STARsolo`
- Outputs: `Aligned.sortedByCoord.out.bam`, `Solo.out`, `counts_unfiltered`

#### STARsolo Barcode Handling
- Removed the linker sequence in R1 so that the `1MM_multi` barcode correction in `STARsolo` can be used
- Barcode & UMI parameters for [`STAR`](https://github.com/alexdobin/STAR/blob/master/docs/STARsolo.md):
```bash
--soloType CB_UMI_Simple \
--soloUMIstart 14 \
--soloUMIlen 7 \
--soloCBstart 1 \
--soloCBlen 14
```

### (u) uTAR pipeline
*U*nnanotated *T*ranscriptionally *A*ctive *R*egions (uTARs)
- Identify unannotated transcriptionally active regions (TARs) using `uTAR_lite`.
- Uses `groHMM` to identify peaks in the STAR alignment outputs, then generates a count matrix.
- ported from https://github.com/mckellardw/uTAR_lite
- Link to paper: https://www.nature.com/articles/s41467-021-22496-3

#### Outputs
- RefFlat file for annotations.
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

