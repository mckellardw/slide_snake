# Borrowed portions of this code from `sockeye` - https://github.com/nanoporetech/sockeye/tree/master


# Merge all the input files into a .fastq
rule ont_1a_merge_formats:
    output:
        MERGED_FQ=temp("{OUTDIR}/{SAMPLE}/ont/tmp/merged.fq.gz"),
    params:
        TMPDIR="{OUTDIR}/{SAMPLE}/ont/tmp",
        ONT_reads=lambda wildcards: ONT[wildcards.SAMPLE],
        CHUNK_SIZE=50,
        OUTPUT_FORMAT="fastq",
    log:
        log="{OUTDIR}/{SAMPLE}/ont/misc_logs/1a_merged.log",
    threads: config["CORES"]
    shell:
        """
        bash scripts/bash/merge_formats_ONT.sh \
            -d "{params.TMPDIR}" \
            -r "{params.ONT_reads}" \
            -c "{params.CHUNK_SIZE}" \
            -o "{output.MERGED_FQ}" \
            -l "{log.log}" \
            -t "{threads}"
        """


# TODO
# Remove super short reads that likely do not have a barcode...
# rule ont_1a_length_filter:
#     input:
#         MERGED_FQ="{OUTDIR}/{SAMPLE}/ont/tmp/merged.fq.gz",
#     input:
#         MERGED_FQ="{OUTDIR}/{SAMPLE}/ont/tmp/merged_filtered.fq.gz",
#     params:
#         MIN_LENGTH=50,
#     threads: config["CORES"]
#     shell:
#         """
#         """

# Remove super short reads that likely do not have a barcode...
# rule ont_1a_pychopper:
#     input:
#         MERGED_FQ="{OUTDIR}/{SAMPLE}/ont/tmp/merged.fq.gz",
#     input:
#         MERGED_FQ="{OUTDIR}/{SAMPLE}/ont/tmp/merged_chopped.fq.gz",
#     params:
#         MIN_LENGTH=100,
#     log:
#         log="{OUTDIR}/{SAMPLE}/ont/misc_logs/1a_pychopper.log",
#     conda:
#         f"{workflow.basedir}/envs/pychopper.yml"
#     resources:
#         mem="16G",
#     threads: 56
#     shell:
#         """
#         """


# borrowed/modified from sockeye (https://github.com/jang1563/sockeye - original ONT github deleted!)
rule ont_1a_call_adapter_scan_v2:
    input:
        FQ="{OUTDIR}/{SAMPLE}/ont/tmp/merged.fq.gz",
    output:
        TSV="{OUTDIR}/{SAMPLE}/ont/adapter_scan.tsv",
        FQ="{OUTDIR}/{SAMPLE}/ont/tmp/merged_stranded.fq.gz",
        # ADAPTERS="{OUTDIR}/{SAMPLE}/ont/adapter_seqs.fasta",
    params:
        BATCH_SIZE=100000,
        ADAPTER1_SEQ="CTACACGACGCTCTTCCGATCT",  #TXG/Curio
        # ADAPTER1_SEQ=lambda w: get_recipe_info(w, "fwd_primer"),
        ADAPTER2_SEQ="ATGTACTCTGCGTTGATACCACTGCTT",  #TXG/Curio
        # ADAPTER2_SEQ="GAGAGAGGAAAGAGAGAGAGAGGG",  #uMRT
        # ADAPTER2_SEQ=lambda w: get_recipe_info(w, "rev_primer"),
        VSEARCH_MIN_ADAPTER_ID=0.7,
    log:
        log="{OUTDIR}/{SAMPLE}/ont/misc_logs/1a_adapter_scan.log",
        err="{OUTDIR}/{SAMPLE}/ont/misc_logs/1a_adapter_scan.err",
    conda:
        f"{workflow.basedir}/envs/adapter_scan.yml"
    resources:
        mem="16G",
    threads: 56
    shell:
        """
        python scripts/py/adapter_scan_vsearch_v2.py \
            --fq_in "{input.FQ}" \
            --fq_out "{output.FQ}" \
            --output_tsv "{output.TSV}" \
            --threads {threads} \
            --batch_size {params.BATCH_SIZE} \
            --adapter1_seq "{params.ADAPTER1_SEQ}" \
            --adapter2_seq "{params.ADAPTER2_SEQ}" \
            --min_adapter_id {params.VSEARCH_MIN_ADAPTER_ID} \
        1> {log.log} \
        2> {log.err}
        """
        # --adapters_fasta "{output.ADAPTERS}" \


# Write lists of read IDs for each adapter type
rule ont_1a_readIDs_by_adapter_type:
    input:
        TSV="{OUTDIR}/{SAMPLE}/ont/adapter_scan.tsv",
        FQ="{OUTDIR}/{SAMPLE}/ont/tmp/merged_stranded.fq.gz",
    output:
        FULL_LEN="{OUTDIR}/{SAMPLE}/ont/adapter_scan_readids/full_len.txt",  # keep
        SINGLE_ADAPTER1="{OUTDIR}/{SAMPLE}/ont/adapter_scan_readids/single_adapter1.txt",  # keep
        # DOUBLE_ADAPTER1="{OUTDIR}/{SAMPLE}/ont/adapter_scan_readids/double_adapter1.txt",  # toss
        # DOUBLE_ADAPTER2="{OUTDIR}/{SAMPLE}/ont/adapter_scan_readids/double_adapter2.txt",  # toss
        # NO_ADAPTERS="{OUTDIR}/{SAMPLE}/ont/adapter_scan_readids/no_adapters.txt",  # toss
        # OTHER="{OUTDIR}/{SAMPLE}/ont/adapter_scan_readids/other.txt",  # toss
        # SINGLE_ADAPTER2="{OUTDIR}/{SAMPLE}/ont/adapter_scan_readids/single_adapter2.txt",  # toss
    resources:
        mem="16G",
    threads: 1
    shell:
        """
        python scripts/py/adapterscan_write_read_id_lists.py \
            --tsv_file_path {input.TSV} \
            --output_directory $(dirname {output.FULL_LEN})
        """

# Summarize adapter_scan results
rule ont_1a_adapter_scan_summary:
    input:
        TSV="{OUTDIR}/{SAMPLE}/ont/adapter_scan.tsv",
    output:
        CSV="{OUTDIR}/{SAMPLE}/ont/misc_logs/1a_adapter_scan_summary.csv",
        PLOT="{OUTDIR}/{SAMPLE}/ont/misc_logs/1a_adapter_scan_summary.pdf",
    params:
        DEVICE="pdf"
    log:
        log="{OUTDIR}/{SAMPLE}/ont/misc_logs/1a_adapter_scan_summary.log",
        err="{OUTDIR}/{SAMPLE}/ont/misc_logs/1a_adapter_scan_summary.err",
    resources:
        mem="8G",
    threads: 1
    conda:
        f"{workflow.basedir}/envs/ggplot2.yml"
    shell:
        """
        Rscript scripts/R/adapter_scan_summary.R \
            -i {input.TSV} \
            -s {output.CSV} \
            -p {output.PLOT} \
            -d {params.DEVICE} \
        1> {log.log} \
        2> {log.err}
        """


# merge lists of useable reads
## FULL_LEN = R1 sequence & TSO sequence
## SINGLE_ADAPTER1 = just R1 sequence - incompletely sequenced
rule ont_1a_merge_scan_lists:
    input:
        FULL_LEN="{OUTDIR}/{SAMPLE}/ont/adapter_scan_readids/full_len.txt",
        SINGLE_ADAPTER1="{OUTDIR}/{SAMPLE}/ont/adapter_scan_readids/single_adapter1.txt",
    output:
        LST="{OUTDIR}/{SAMPLE}/ont/adapter_scan_readids/keep.txt",
    resources:
        mem="8G",
    threads: 1
    shell:
        """
        cat {input.FULL_LEN} {input.SINGLE_ADAPTER1} > {output.LST}
        """


# TODO- add more functionality for other read/adapter types to salvage imperfect reads
rule ont_1a_subset_fastq_by_adapter_type:
    input:
        FQ="{OUTDIR}/{SAMPLE}/ont/tmp/merged_stranded.fq.gz",
        LST="{OUTDIR}/{SAMPLE}/ont/adapter_scan_readids/keep.txt",
        # FULL_LEN = "{OUTDIR}/{SAMPLE}/ont/adapter_scan_readids/full_len.txt",
        # SINGLE_ADAPTER1 = "{OUTDIR}/{SAMPLE}/ont/adapter_scan_readids/single_adapter1.txt",
    output:
        FQ=temp("{OUTDIR}/{SAMPLE}/ont/tmp/merged_stranded.fq"),
        # FULL_LEN = "{OUTDIR}/{SAMPLE}/ont/tmp/adapter_scan_readids/full_len.fq.gz",
        # SINGLE_ADAPTER1 = "{OUTDIR}/{SAMPLE}/ont/tmp/adapter_scan_readids/single_adapter1.fq.gz",
        FQ_ADAPTER="{OUTDIR}/{SAMPLE}/ont/tmp/merged_adapter.fq",
    resources:
        mem="16G",
    threads: 1
    log:
        err="{OUTDIR}/{SAMPLE}/ont/misc_logs/1a_subset_fastq_by_adapter_type.err",
    conda:
        f"{workflow.basedir}/envs/seqkit.yml"
    shell:
        """
        mkdir -p $(dirname {output.FQ})            
        zcat {input.FQ} > {output.FQ} 

        seqkit grep \
            -f {input.LST} \
            {output.FQ} \
        > {output.FQ_ADAPTER} \
        2> {log.err}
        """


# Compress the merged fastq
rule ont_1a_compress_merged_fq:
    input:
        FQ="{OUTDIR}/{SAMPLE}/ont/tmp/merged_adapter.fq",
    output:
        FQ="{OUTDIR}/{SAMPLE}/ont/tmp/merged_adapter.fq.gz",
    resources:
        mem="8G",
    threads: config["CORES"]
    shell:
        """
        pigz -p{threads} {input.FQ}
        """


# Split reads in the poly(T) stretch and rev-comp the "R1" sequence
##TODO: add read length bounds for R1 based on barcode construct to reduce incorrect split sites across reads
rule ont_1a_split_fastq_to_R1_R2:
    input:
        FQ="{OUTDIR}/{SAMPLE}/ont/tmp/merged_adapter.fq.gz",
    output:
        R1_FQ="{OUTDIR}/{SAMPLE}/ont/tmp/merged_adapter_R1.fq.gz",
        R2_FQ="{OUTDIR}/{SAMPLE}/ont/tmp/merged_adapter_R2.fq.gz",
        AMBIG_FQ="{OUTDIR}/{SAMPLE}/ont/tmp/merged_adapter_ambiguous.fq.gz",
    params:
        # ANCHOR_SEQ=lambda w: get_recipe_info(w, "fwd_primer"),
        # ANCHOR_SEQ="CTACACGACGCTCTTCCGATCT",  #TXG/Curio
        ANCHOR_SEQ="ACGCTCTTCCGATCT",  #TXG/Curio
        SPLIT_SEQ="T" * 8,
        SPLIT_OFFSET=8,  # offset from 3' end of split seq on which to split
        MAX_OFFSET=200,
    resources:
        mem="16G",
    # threads:  config["CORES"]
    log:
        log="{OUTDIR}/{SAMPLE}/ont/misc_logs/1a_read_split.log",
        err="{OUTDIR}/{SAMPLE}/ont/misc_logs/1a_read_split.err",
    conda:
        f"{workflow.basedir}/envs/parasail.yml"
    shell:
        """
        python scripts/py/fastq_split_reads_parallelized_v3.py --fq_in {input.FQ} \
            --anchor_seq {params.ANCHOR_SEQ} \
            --split_seq {params.SPLIT_SEQ} \
            --split_offset {params.SPLIT_OFFSET} \
            --max_offset {params.MAX_OFFSET} \
            --threads {threads} \
        1> {log.log} \
        2> {log.err}
        """
