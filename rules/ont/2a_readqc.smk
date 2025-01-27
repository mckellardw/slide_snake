# fastqc before adapter scan
rule ont_2a_readQC_0_rawInput:
    input:
        FQ="{OUTDIR}/{SAMPLE}/ont/tmp/merged.fq.gz",
    output:
        TSV="{OUTDIR}/{SAMPLE}/ont/readqc/0_rawInput/merged_qc.tsv",
    params:
        CHUNK_SIZE=100000,  # number of reads to handle in each chunk
    log:
        log="{OUTDIR}/{SAMPLE}/ont/readqc/0_rawInput/merged_qc.log",
    resources:
        mem="8G",
    threads: config["CORES"]
    shell:
        """
        python scripts/py/fastq_readqc.py \
            {input.FQ} \
            {output.TSV} \
            --cores {threads} \
            --chunk_size {params.CHUNK_SIZE} \
        > {log.log}
        """


# fastqc before trimming
rule ont_2a_readQC_1_preCutadapt:
    input:
        FQ="{OUTDIR}/{SAMPLE}/ont/tmp/merged_adapter_{READ}.fq.gz",
    output:
        TSV="{OUTDIR}/{SAMPLE}/ont/readqc/1_preCutadapt/{READ}_qc.tsv",
    params:
        CHUNK_SIZE=500000,  # number of reads to handle in each chunk
    log:
        log="{OUTDIR}/{SAMPLE}/ont/readqc/1_preCutadapt/{READ}_qc.log",
    resources:
        mem="8G",
    threads: config["CORES"]
    shell:
        """
        python scripts/py/fastq_readqc.py \
            {input.FQ} \
            {output.TSV} \
            --cores {threads} \
            --chunk_size {params.CHUNK_SIZE} \
        > {log.log}
        """


# fastqc after cutadapt trimming
rule ont_2a_readQC_2_postCutadapt:
    input:
        FQ="{OUTDIR}/{SAMPLE}/ont/tmp/cut_{READ}.fq.gz",
    output:
        TSV="{OUTDIR}/{SAMPLE}/ont/readqc/2_postCutadapt/{READ}_qc.tsv",
    params:
        CHUNK_SIZE=500000,
    log:
        log="{OUTDIR}/{SAMPLE}/ont/readqc/2_postCutadapt/{READ}_qc.log",
    resources:
        mem="8G",
    threads: config["CORES"]
    shell:
        """
        python scripts/py/fastq_readqc.py \
            {input.FQ} \
            {output.TSV} \
            --cores {threads} \
            --chunk_size {params.CHUNK_SIZE} \
        > {log.log}
        """


# ┌────┬──────┬───────────────────────────────────────────────────────┐
# │Tag │ Type │                      Description                      │
# ├────┼──────┼───────────────────────────────────────────────────────┤
# │ tp │  A   │ Type of aln: P/primary, S/secondary and I,i/inversion │
# │ cm │  i   │ Number of minimizers on the chain                     │
# │ s1 │  i   │ Chaining score                                        │
# │ s2 │  i   │ Chaining score of the best secondary chain            │
# │ NM │  i   │ Total number of mismatches and gaps in the alignment  │
# │ MD │  Z   │ To generate the ref sequence in the alignment         │
# │ AS │  i   │ DP alignment score                                    │
# │ ms │  i   │ DP score of the max scoring segment in the alignment  │
# │ nn │  i   │ Number of ambiguous bases in the alignment            │
# │ ts │  A   │ Transcript strand (splice mode only)                  │
# │ cg │  Z   │ CIGAR string (only in PAF)                            │
# │ cs │  Z   │ Difference string                                     │
# └────┴──────┴───────────────────────────────────────────────────────┘
rule ont_2a_readQC_3_bam:
    input:
        BAM="{OUTDIR}/{SAMPLE}/ont/minimap2/{RECIPE}/sorted_gn_cb_ub.bam",
        BAI="{OUTDIR}/{SAMPLE}/ont/minimap2/{RECIPE}/sorted_gn_cb_ub.bam.bai",
    output:
        TSV="{OUTDIR}/{SAMPLE}/ont/readqc/3_aligned/{RECIPE}_qc.tsv",
    params:
        TAGS="AS NM GN CB UR",
        CHUNK_SIZE=500000,
    log:
        log="{OUTDIR}/{SAMPLE}/ont/readqc/3_aligned/{RECIPE}_qc.log",
        err="{OUTDIR}/{SAMPLE}/ont/readqc/3_aligned/{RECIPE}_qc.err",
    resources:
        mem="8G",
    threads: 1
    # threads=config["CORES"],
    shell:
        """
        python scripts/py/bam_readqc.py \
            --tags {params.TAGS} \
            --chunk-size {params.CHUNK_SIZE} \
            --bam_file {input.BAM} \
            --tsv_file {output.TSV} \
        1> {log.log} \
        2> {log.err}
        """


# Grab the first million reads...
rule ont_2a_readQC_downsample:
    input:
        TSV="{OUTDIR}/{SAMPLE}/ont/readqc/{TRIM}/{READ}_qc.tsv",
    output:
        TSV="{OUTDIR}/{SAMPLE}/ont/readqc/{TRIM}/{READ}_qc_downSampled.tsv",
    params:
        N_READS=config["N_READS_READQC"],
    resources:
        mem="4G",
    threads: 1
    shell:
        """
        head -n 1 {input.TSV} > {output.TSV} && tail -n +2 {input.TSV} | shuf -n {params.N_READS} >> {output.TSV}
        """


# Summary plot
rule ont_2a_readQC_summaryplot:
    input:
        TSV="{OUTDIR}/{SAMPLE}/ont/readqc/{TRIM}/{READ}_qc_downSampled.tsv",
    output:
        IMG="{OUTDIR}/{SAMPLE}/ont/readqc/{TRIM}/{READ}_qc.png",
    resources:
        mem="8G",
    threads: 1
    log:
        log="{OUTDIR}/{SAMPLE}/short_read/readqc/{TRIM}/{READ}_qc.log",
        err="{OUTDIR}/{SAMPLE}/short_read/readqc/{TRIM}/{READ}_qc.err",
    conda:
        f"{workflow.basedir}/envs/ggplot2.yml"
    shell:
        """
        Rscript scripts/R/readqc_summary.R \
            -f {input.TSV} \
            -o {output.IMG} \
        1> {log.log} \
        2> {log.err}
        """


# Compress TSV files using pigz
rule ont_2a_readQC_compress:
    input:
        TSV="{OUTDIR}/{SAMPLE}/ont/readqc/{TRIM}/{READ}_qc.tsv",
    output:
        GZ="{OUTDIR}/{SAMPLE}/ont/readqc/{TRIM}/{READ}_qc.tsv.gz",
    resources:
        mem="4G",
    threads: 1
    shell:
        """
        pigz -c {input.TSV} > {output.GZ}
        """
