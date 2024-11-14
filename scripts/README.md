# scripts

## Python Scripts

- `py/adapter_scan_vsearch_v2.py` - Scans for adapter sequences in reads using VSEARCH and processes the results.
- `py/adapterscan_write_read_id_lists.py` - Processes a TSV file to create text files for each unique "lab" value containing corresponding "read_id" values.
- `py/add_sequence_as_tag.py` - Adds sequences from a FASTQ file as tags in a BAM file.
- `py/bam_readID2tags.py` - Parses a BAM file, extracts the read ID, and moves the last two strings in the ID to tags.
- `py/bam_readqc.py` - Calculates various quality control metrics for reads in a BAM file.
- `py/cache_umitools_to_h5ad.py` - Loads a UMI-tools count matrix and spatial data, then saves it as an AnnData object.
- `py/fastq_call_bc_umi_from_adapter_v2.py` - Calls barcodes and UMIs from adapter sequences in FASTQ files.
- `py/fastq_internal_adapter_trim_R1_v2.py` - Removes internal adapter sequences from R1 reads in FASTQ files.
- `py/fastq_readqc.py` - Calculates various quality control metrics for reads in a FASTQ file.
- `py/fastqs2bam.py` - Converts paired FASTQ files to an unaligned BAM file, adding barcode information from Read 1 as a tag and quality strings as tags.
- `py/long2mtx.py` - Converts a TSV file to a sparse matrix in mtx format.
- `py/tsv2tag.py` - Adds gene assignment from a TSV file to a BAM file as a tag.
- `py/tsv_bc_correction_parallelized.py` - Processes a TSV file to correct barcodes based on a whitelist, in parallel.
- `py/tsv_bc_correction_summary.py` - Calculates statistics for barcode correction results and writes them to a file.

## Shell Scripts

- `bam_dedup.sh` - Deduplicates BAM files aligned with STARsolo, processing chromosome by chromosome.
- `bam_dedupByChr.sh` - Deduplicates BAM files aligned with STARsolo, processing chromosome by chromosome with additional options.
- `kb.sh` - Runs kallisto/bustools for RNA-seq data processing.
- `splitNfqs.py` - Splits a FASTQ file into N chunks using a parallelized `sed` command.

## R Scripts

- `adapter_scan_summary.R` - Summarizes adapter scan results and generates plots.
- `cutadapt_summary.R` - Summarizes Cutadapt results and generates plots.
- `readqc_summary.R` - Summarizes read quality control results and generates plots.