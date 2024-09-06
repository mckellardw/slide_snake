# Filter out rRNA reads w/ bwa alignment
# VASAseq implementation - https://github.com/anna-alemany/VASAseq/blob/main/mapping/ribo-bwamem.sh
# Align to rRNA ref w/ `bwa mem` for cleaner/faster rRNA filtering
##TODO incorporate VASAseq style "long"/short read handling with multiple align steps


# TODO- refactor to incorporate internal trimming options into rRNA filtering
rule ilmn_2a_bwa_rRNA_align:
    input:
        R2_FQ="{OUTDIR}/{SAMPLE}/tmp/twiceCut_R2.fq.gz",
        BC_WHITELIST="{OUTDIR}/{SAMPLE}/bc/whitelist.txt",
        BC_1="{OUTDIR}/{SAMPLE}/bc/whitelist_1.txt",
        BC_2="{OUTDIR}/{SAMPLE}/bc/whitelist_2.txt",
        BC_ADAPTER="{OUTDIR}/{SAMPLE}/bc/whitelist_adapter.txt",
    output:
        BAM1=temp("{OUTDIR}/{SAMPLE}/rRNA/bwa/aligned.bam"),
        BAM2="{OUTDIR}/{SAMPLE}/rRNA/bwa/aligned_sorted.bam",
        # R1_FQ_BWA_FILTERED  = '{OUTDIR}/{SAMPLE}/rRNA/bwa/noRibo_R1.fq',
        R2_FQ_BWA_FILTERED="{OUTDIR}/{SAMPLE}/rRNA/bwa/noRibo_R2.fq",
    params:
        # MEMLIMIT = config['MEMLIMIT'],
        BWA_REF=lambda w: get_bwa_ref(w, mode="rRNA"),
        MIN_ALIGNSCORE=40,
    log:
        log="{OUTDIR}/{SAMPLE}/rRNA/bwa/bwa_mem.log",
    resources:
        mem="96G",
        time="2:00:00",
    threads: config["CORES"]
    conda:
        f"{workflow.basedir}/envs/bwa.yml"
    shell:
        """
        mkdir -p $(dirname {output.BAM1})

        bwa-mem2 mem \
            -t {threads} \
            {params.BWA_REF} \
            {input.R2_FQ} \
        1> {output.BAM1} \
        2> {log.log} \
        
        samtools sort \
            -@ {threads} \
            -O BAM \
            {output.BAM1} \
        > {output.BAM2} 

        samtools view \
            -h {output.BAM2} \
        | awk \
            -v quality={params.MIN_ALIGNSCORE} \
            -f scripts/awk/bam_lowPassMAPQFilter.awk \
        | samtools fastq \
            {output.BAM2} \
        > {output.R2_FQ_BWA_FILTERED} 
        """


rule ilmn_2a_bwa_rRNA_filter_R1:
    input:
        R1_FQ="{OUTDIR}/{SAMPLE}/tmp/twiceCut_R1.fq.gz",
        R2_FQ="{OUTDIR}/{SAMPLE}/tmp/twiceCut_R2.fq.gz",
        R2_FQ_BWA_FILTERED="{OUTDIR}/{SAMPLE}/rRNA/bwa/noRibo_R2.fq",
    output:
        R1_FQ=temp("{OUTDIR}/{SAMPLE}/tmp/cut_R1.fq"),
        R1_FQ_BWA_FILTERED="{OUTDIR}/{SAMPLE}/rRNA/bwa/noRibo_R1.fq",
        rRNA_LIST="{OUTDIR}/{SAMPLE}/rRNA/bwa/rRNA_readID.list",
    resources:
        mem="64G",
    threads: config["CORES"]
    shell:
        """
        cat {input.R2_FQ_BWA_FILTERED} | awk -f scripts/awk/fq_readHeader.awk - > {output.rRNA_LIST}

        zcat {input.R1_FQ} > {output.R1_FQ} 
        
        seqtk subseq {output.R1_FQ} {output.rRNA_LIST} > {output.R1_FQ_BWA_FILTERED}
        """


rule ilmn_2a_bwa_rRNA_compress_unmapped:
    input:
        R1_FQ_BWA_FILTERED="{OUTDIR}/{SAMPLE}/rRNA/bwa/noRibo_R1.fq",
        R2_FQ_BWA_FILTERED="{OUTDIR}/{SAMPLE}/rRNA/bwa/noRibo_R2.fq",
    output:
        R1_FQ_BWA_FILTERED="{OUTDIR}/{SAMPLE}/rRNA/bwa/noRibo_R1.fq.gz",
        R2_FQ_BWA_FILTERED="{OUTDIR}/{SAMPLE}/rRNA/bwa/noRibo_R2.fq.gz",
    resources:
        mem="16G",
    threads: config["CORES"]
    shell:
        """
        pigz -p{threads} {input}
        """


#  Run fastqc on unmapped reads;
rule ilmn_2a_bwa_rRNA_filtered_fastqc:
    input:
        FQ="{OUTDIR}/{SAMPLE}/rRNA/bwa/noRibo_{READ}.fq.gz",
    output:
        FQC_DIR=directory("{OUTDIR}/{SAMPLE}/fastqc/rRNA_bwa_{READ}"),
    params:
        adapters=config["FASTQC_ADAPTERS"],
    resources:
        mem="16G",
    threads: config["CORES"]
    conda:
        f"{workflow.basedir}/envs/fastqc.yml"
    shell:
        """
        mkdir -p {output.FQC_DIR}

        fastqc \
            -o {output.FQC_DIR} \
            -t {threads} \
            -a {params.adapters} \
            {input.FQ}
        """
