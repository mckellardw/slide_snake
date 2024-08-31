# Borrowed portions of this code from `sockeye` - https://github.com/nanoporetech/sockeye/tree/master


# TODO- rewrite as a python script...
rule ont_1a_merge_formats:
    output:
        MERGED_FQ=temp("{OUTDIR}/{SAMPLE}/tmp/ont/merged.fq.gz"),
    params:
        TMPDIR="{OUTDIR}/{SAMPLE}/tmp/ont",
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


## Deprecated
# rule ont_call_adapter_scan:
#     input:
#         FQ="{OUTDIR}/{SAMPLE}/tmp/ont/merged.fq.gz",
#     output:
#         TSV="{OUTDIR}/{SAMPLE}/ont/adapter_scan.tsv",
#         FQ="{OUTDIR}/{SAMPLE}/tmp/ont/merged_stranded.fq.gz",
#     params:
#         # batch_size = config["READ_STRUCTURE_BATCH_SIZE"],
#         # KIT="3prime",  #['3prime', '5prime', 'multiome']
#         KIT="uMRT",
#     log:
#         log="{OUTDIR}/{SAMPLE}/ont/misc_logs/adapter_scan.log",
#     conda:
#         f"{workflow.basedir}/envs/adapter_scan.yml"
#     resources:
#         mem="16G",
#     threads: 56
#     shell:
#         """
#         python scripts/py/adapter_scan_vsearch.py \
#             --kit {params.KIT} \
#             --output_fastq {output.FQ} \
#             --output_tsv {output.TSV} \
#             -t {threads} \
#             {input.FQ} \
#         2>&1 | tee {log.log}
#         """


# borrowed/modified from sockeye (https://github.com/jang1563/sockeye - original ONT github deleted!)
rule ont_1a_call_adapter_scan_v2:
    input:
        FQ="{OUTDIR}/{SAMPLE}/tmp/ont/merged.fq.gz",
    output:
        TSV="{OUTDIR}/{SAMPLE}/ont/adapter_scan.tsv",
        FQ="{OUTDIR}/{SAMPLE}/tmp/ont/merged_stranded.fq.gz",
        # ADAPTERS="{OUTDIR}/{SAMPLE}/ont/adapter_seqs.fasta",
    params:
        BATCH_SIZE=100000,
        ADAPTER1_SEQ="CTACACGACGCTCTTCCGATCT",  #TXG/Curio
        # ADAPTER2_SEQ="ATGTACTCTGCGTTGATACCACTGCTT", #TXG/Curio
        ADAPTER2_SEQ="GAGAGAGGAAAGAGAGAGAGAGGG",  #uMRT
        VSEARCH_MIN_ADAPTER_ID=0.7,
    log:
        log="{OUTDIR}/{SAMPLE}/ont/misc_logs/1a_adapter_scan.log",
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
        2>&1 | tee {log.log}
        """
        # --adapters_fasta "{output.ADAPTERS}" \


# Write lists of read IDs for each adapter type
rule ont_1a_readIDs_by_adapter_type:
    input:
        TSV="{OUTDIR}/{SAMPLE}/ont/adapter_scan.tsv",
        FQ="{OUTDIR}/{SAMPLE}/tmp/ont/merged_stranded.fq.gz",
    output:
        FULL_LEN="{OUTDIR}/{SAMPLE}/ont/adapter_scan_readids/full_len.txt",  # keep
        SINGLE_ADAPTER1="{OUTDIR}/{SAMPLE}/ont/adapter_scan_readids/single_adapter1.txt",  # keep
        DOUBLE_ADAPTER1="{OUTDIR}/{SAMPLE}/ont/adapter_scan_readids/double_adapter1.txt",  # toss
        DOUBLE_ADAPTER2="{OUTDIR}/{SAMPLE}/ont/adapter_scan_readids/double_adapter2.txt",  # toss
        NO_ADAPTERS="{OUTDIR}/{SAMPLE}/ont/adapter_scan_readids/no_adapters.txt",  # toss
        OTHER="{OUTDIR}/{SAMPLE}/ont/adapter_scan_readids/other.txt",  # toss
        SINGLE_ADAPTER2="{OUTDIR}/{SAMPLE}/ont/adapter_scan_readids/single_adapter2.txt",  # toss
    resources:
        mem="16G",
    threads: config["CORES"]
    shell:
        """
        python scripts/py/adapterscan_write_read_id_lists.py \
            --tsv_file_path {input.TSV} \
            --output_directory $(dirname {output.FULL_LEN})
        """


# Write lists of read IDs for each adapter type
rule ont_1a_adapter_scan_results:
    input:
        FULL_LEN="{OUTDIR}/{SAMPLE}/ont/adapter_scan_readids/full_len.txt",
        SINGLE_ADAPTER1="{OUTDIR}/{SAMPLE}/ont/adapter_scan_readids/single_adapter1.txt",
    output:
        LOG="{OUTDIR}/{SAMPLE}/ont/misc_logs/1a_adapter_scan_results.txt",
    resources:
        mem="8G",
    threads: 1
    shell:
        """
        dir_path=$(dirname {input.FULL_LEN})

        for file in "$dir_path"/*.txt; do
            echo "$(basename $file)"\t"$(wc -l <"$file")" >> {output.LOG}
        done
        """


# merge lists of useable reads
## FULL_LEN = R1 sequence & TSO sequence
## SINGLE_ADAPTER1 = just R1 sequence - incompletely sequenced
rule ont_1a_merge_scan_lists:
    input:
        LOG="{OUTDIR}/{SAMPLE}/ont/misc_logs/1a_adapter_scan_results.txt",
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
        FQ="{OUTDIR}/{SAMPLE}/tmp/ont/merged_stranded.fq.gz",
        LST="{OUTDIR}/{SAMPLE}/ont/adapter_scan_readids/keep.txt",
        # FULL_LEN = "{OUTDIR}/{SAMPLE}/ont/adapter_scan_readids/full_len.txt",
        # SINGLE_ADAPTER1 = "{OUTDIR}/{SAMPLE}/ont/adapter_scan_readids/single_adapter1.txt",
    output:
        FQ=temp("{OUTDIR}/{SAMPLE}/tmp/ont/merged_stranded.fq"),
        # FULL_LEN = "{OUTDIR}/{SAMPLE}/tmp/ont/adapter_scan_readids/full_len.fq.gz",
        # SINGLE_ADAPTER1 = "{OUTDIR}/{SAMPLE}/tmp/ont/adapter_scan_readids/single_adapter1.fq.gz",
        FQ_ADAPTER="{OUTDIR}/{SAMPLE}/tmp/ont/adapter_scan_readids/merged_adapter.fq",
    resources:
        mem="16G",
    threads: 1
    log:
        log="{OUTDIR}/{SAMPLE}/ont/misc_logs/1a_subset_fastq_by_adapter_type.log",
    shell:
        """
        mkdir -p $(dirname {output.FQ})            
        zcat {input.FQ} > {output.FQ} 

        seqtk subseq \
            {output.FQ} \
            {input.LST} \
        > {output.FQ_ADAPTER} \
        2> {log.log}
        """


rule ont_1a_compress_merged_fq:
    input:
        FQ="{OUTDIR}/{SAMPLE}/tmp/ont/adapter_scan_readids/merged_adapter.fq",
    output:
        FQ="{OUTDIR}/{SAMPLE}/tmp/ont/adapter_scan_readids/merged_adapter.fq.gz",
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
        FQ="{OUTDIR}/{SAMPLE}/tmp/ont/adapter_scan_readids/merged_adapter.fq.gz",
    output:
        R1_FQ="{OUTDIR}/{SAMPLE}/tmp/ont/adapter_scan_readids/merged_adapter_R1.fq.gz",
        R2_FQ="{OUTDIR}/{SAMPLE}/tmp/ont/adapter_scan_readids/merged_adapter_R2.fq.gz",
    params:
        ADAPTER="T" * 8,
    resources:
        mem="16G",
    threads: config["CORES"]
    log:
        log="{OUTDIR}/{SAMPLE}/ont/misc_logs/1a_read_split.log",
    run:
        # for ADAPTER in input.ADAPTER_TYPES: #TODO- broaden to other read types, beyond full_len
        shell(
            f"""
            python scripts/py/fastq_split_reads_parallelized.py --fq_in {input.FQ} \
                --split_seq {params.ADAPTER} \
                --threads {threads} \
            2> {log.log}
            """
        )
        # --min_R1_length {} \
        # --max_R1_length {} \
