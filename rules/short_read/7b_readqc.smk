# fastqc before adapter scan
rule ilmn_7b_readQC_0_rawInput:
    input:
        FQ="{OUTDIR}/{SAMPLE}/short_read/tmp/merged_{READ}.fq.gz",
    output:
        TSV="{OUTDIR}/{SAMPLE}/short_read/readqc/0_rawInput/{READ}_qc.tsv",
    params:
        CHUNK_SIZE=100000,  # number of reads to handle in each chunk
    log:
        log="{OUTDIR}/{SAMPLE}/short_read/readqc/0_rawInput/{READ}_qc.log",
        err="{OUTDIR}/{SAMPLE}/short_read/readqc/0_rawInput/{READ}_qc.err",
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
        > {log.log} \
        2> {log.err}
        """


# fastqc before trimming
rule ilmn_7b_readQC_1_preCutadapt:
    input:
        FQ="{OUTDIR}/{SAMPLE}/short_read/tmp/merged_{READ}.fq.gz",
    output:
        TSV="{OUTDIR}/{SAMPLE}/short_read/readqc/1_preCutadapt/{READ}_qc.tsv",
    params:
        CHUNK_SIZE=500000,  # number of reads to handle in each chunk
    log:
        log="{OUTDIR}/{SAMPLE}/short_read/readqc/1_preCutadapt/{READ}_qc.log",
        err="{OUTDIR}/{SAMPLE}/short_read/readqc/1_preCutadapt/{READ}_qc.err",
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
        > {log.log} \
        2> {log.err}
        """


# fastqc after cutadapt trimming
rule ilmn_7b_readQC_2_postCutadapt:
    input:
        FQ="{OUTDIR}/{SAMPLE}/short_read/tmp/cut_{READ}.fq.gz",
    output:
        TSV="{OUTDIR}/{SAMPLE}/short_read/readqc/2_postCutadapt/{READ}_qc.tsv",
    params:
        CHUNK_SIZE=500000,
    log:
        log="{OUTDIR}/{SAMPLE}/short_read/readqc/2_postCutadapt/{READ}_qc.log",
        err="{OUTDIR}/{SAMPLE}/short_read/readqc/2_postCutadapt/{READ}_qc.err",
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
        > {log.log} \
        2> {log.err}
        """


# fastqc after second cutadapt trimming
rule ilmn_7b_readQC_3_twiceCutadapt:
    input:
        FQ="{OUTDIR}/{SAMPLE}/short_read/tmp/twiceCut_{READ}.fq.gz",
    output:
        TSV="{OUTDIR}/{SAMPLE}/short_read/readqc/3_twiceCutadapt/{READ}_qc.tsv",
    params:
        CHUNK_SIZE=500000,
    log:
        log="{OUTDIR}/{SAMPLE}/short_read/readqc/3_twiceCutadapt/{READ}_qc.log",
        err="{OUTDIR}/{SAMPLE}/short_read/readqc/3_twiceCutadapt/{READ}_qc.err",
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
        > {log.log} \
        2> {log.err}
        """


# https://github.com/alexdobin/STAR/blob/master/docs/STARsolo.md
# CR/UR: raw (uncorrected) CellBarcode/UMI
# CY/UY: quality score for CellBarcode/UMI
# GX/GN: for gene ID/names
# sS/sQ: for sequence/quality combined CellBarcode and UMI; sM for barcode match status.
# CB/UB: corrected CellBarcode/UMI.
rule ilmn_7b_readQC_3_bam:
    input:
        BAM="{OUTDIR}/{SAMPLE}/short_read/STARsolo/{RECIPE}/Aligned.sortedByCoord.out.bam",
        BAI="{OUTDIR}/{SAMPLE}/short_read/STARsolo/{RECIPE}/Aligned.sortedByCoord.out.bam.bai",
    output:
        TSV="{OUTDIR}/{SAMPLE}/short_read/readqc/4_aligned/{RECIPE}_qc.tsv",
    params:
        TAGS="NH HI nM AS UB CB GX GN sS sQ sM",
        CHUNK_SIZE=500000,
    log:
        log="{OUTDIR}/{SAMPLE}/short_read/readqc/4_aligned/{RECIPE}_qc.log",
        err="{OUTDIR}/{SAMPLE}/short_read/readqc/4_aligned/{RECIPE}_qc.err",
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
rule ilmn_7b_readQC_downsample:
    input:
        TSV="{OUTDIR}/{SAMPLE}/short_read/readqc/{TRIM}/{READ}_qc.tsv",
    output:
        TSV="{OUTDIR}/{SAMPLE}/short_read/readqc/{TRIM}/{READ}_qc_downSampled.tsv",
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
rule ilmn_7b_readQC_summaryplot:
    input:
        TSV="{OUTDIR}/{SAMPLE}/short_read/readqc/{TRIM}/{READ}_qc_downSampled.tsv",
    output:
        IMG="{OUTDIR}/{SAMPLE}/short_read/readqc/{TRIM}/{READ}_qc.png",
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
rule ilmn_7b_readQC_compress:
    input:
        TSV="{OUTDIR}/{SAMPLE}/short_read/readqc/{TRIM}/{READ}_qc.tsv",
    output:
        GZ="{OUTDIR}/{SAMPLE}/short_read/readqc/{TRIM}/{READ}_qc.tsv.gz",
    resources:
        mem="4G",
    threads: 1
    shell:
        """
        pigz -c {input.TSV} > {output.GZ}
        """
