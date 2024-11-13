# Filter out rRNA reads w/ bwa alignment
# VASAseq implementation - https://github.com/anna-alemany/VASAseq/blob/main/mapping/ribo-bwamem.sh
# Align to rRNA ref w/ `bwa mem` for cleaner/faster rRNA filtering
##TODO incorporate VASAseq style "long"/short read handling with multiple align steps
rule ilmn_2a_bwa_rRNA_align:
    input:
        R2_FQ="{OUTDIR}/{SAMPLE}/short_read/tmp/twiceCut_R2.fq.gz",
    output:
        BAM1=temp("{OUTDIR}/{SAMPLE}/short_read/rRNA/bwa/aligned.bam"),
        BAM2="{OUTDIR}/{SAMPLE}/short_read/rRNA/bwa/aligned_sorted.bam",
        R2_FQ_BWA_FILTERED="{OUTDIR}/{SAMPLE}/short_read/rRNA/bwa/noRibo_R2.fq",
    params:
        BWA_REF=lambda w: get_bwa_ref(w, mode="rRNA"),
        MIN_ALIGNSCORE=40,
    log:
        log="{OUTDIR}/{SAMPLE}/short_read/rRNA/bwa/bwa_mem.log",
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
        2> {log.log}

        samtools sort \
            -@ {threads} \
            -O BAM \
            -o {output.BAM2} \
            {output.BAM1}

        samtools view \
            -h {output.BAM2} \
        | awk \
            -v quality={params.MIN_ALIGNSCORE} \
            -f scripts/awk/bam_lowPassMAPQFilter.awk \
        | samtools fastq \
        > {output.R2_FQ_BWA_FILTERED}
        """


# Generate list of read IDs to keep ("no ribo") from the filtered R2 file
rule ilmn_2a_bwa_rRNA_get_noRibo_list:
    input:
        R2_FQ_BWA_FILTERED="{OUTDIR}/{SAMPLE}/short_read/rRNA/bwa/noRibo_R2.fq",
    output:
        rRNA_LIST=temp("{OUTDIR}/{SAMPLE}/short_read/rRNA/bwa/rRNA_readID.list"),
    threads: 1
    shell:
        """
        cat {input.R2_FQ_BWA_FILTERED} | awk -f scripts/awk/fq_readHeader.awk - > {output.rRNA_LIST}
        """

# Filter R1 reads using the generated list of read IDs
rule ilmn_2a_bwa_rRNA_filter_R1:
    input:
        R1_FQ="{OUTDIR}/{SAMPLE}/short_read/tmp/twiceCut_{TRIM}_R1.fq.gz",
        rRNA_LIST="{OUTDIR}/{SAMPLE}/short_read/rRNA/bwa/rRNA_readID.list",
    output:
        R1_FQ_BWA_FILTERED="{OUTDIR}/{SAMPLE}/short_read/rRNA/bwa/noRibo_{TRIM}_R1.fq",
    resources:
        mem="64G",
    threads: config["CORES"]
    conda:
        f"{workflow.basedir}/envs/seqkit.yml"
    shell:
        """
        seqkit grep -f {input.rRNA_LIST} {input.R1_FQ} -o {output.R1_FQ_BWA_FILTERED}
        """


# Compress previous outputs
rule ilmn_2a_bwa_rRNA_compress_unmapped:
    input:
        R1_FQ_BWA_FILTERED="{OUTDIR}/{SAMPLE}/short_read/rRNA/bwa/noRibo_{READ}.fq",
    output:
        R1_FQ_BWA_FILTERED="{OUTDIR}/{SAMPLE}/short_read/rRNA/bwa/noRibo_{READ}.fq.gz",
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
        FQ="{OUTDIR}/{SAMPLE}/short_read/rRNA/bwa/noRibo_{READ}.fq.gz",
    output:
        FQC_DIR=directory("{OUTDIR}/{SAMPLE}/fastqc/rRNA_bwa/{READ}"),
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
