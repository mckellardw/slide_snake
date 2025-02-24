# Filter out rRNA reads w/ bwa alignment
# VASAseq implementation - https://github.com/anna-alemany/VASAseq/blob/main/mapping/ribo-bwamem.sh
##TODO incorporate VASAseq style "long"/short read handling with multiple align steps
# Align to rRNA ref; keep reads with low alignment score (no_rRNA_R2.fq)
# Tag	Description
# AS:i:	Alignment score (higher is better)
# XS:i:	Suboptimal alignment score (for secondary alignments)
# NM:i:	Edit distance (number of mismatches and indels)
# MD:Z:	Mismatch string (indicates mismatches against the reference)
rule ilmn_2a_bwa_rRNA_align:
    input:
        R2_FQ="{OUTDIR}/{SAMPLE}/short_read/tmp/twiceCut_R2.fq.gz",
    output:
        BAM1=temp("{OUTDIR}/{SAMPLE}/short_read/rRNA/bwa/aligned.sam"),
        BAM2="{OUTDIR}/{SAMPLE}/short_read/rRNA/bwa/aligned_sorted.bam",
        R2_FQ_BWA_FILTERED="{OUTDIR}/{SAMPLE}/short_read/rRNA/bwa/no_rRNA_R2.fq",
    params:
        BWA_REF=lambda w: get_bwa_ref(w, mode="rRNA"),
        MIN_ALIGNSCORE=40,
    log:
        log="{OUTDIR}/{SAMPLE}/short_read/rRNA/bwa/bwa_mem.log",
        err="{OUTDIR}/{SAMPLE}/short_read/rRNA/bwa/bwa_mem.err",
    resources:
        mem="96G",
        time="2:00:00",
    threads: config["CORES"]
    conda:
        f"{workflow.basedir}/envs/bwa.yml"
    shell:
        """
        bash scripts/bash/fq_bwa_rRNA_align.sh \
            -r {input.R2_FQ} \
            -a {output.BAM1} \
            -s {output.BAM2} \
            -n {output.R2_FQ_BWA_FILTERED} \
            -f {params.BWA_REF} \
            -q {params.MIN_ALIGNSCORE} \
            -t {threads} \
        1> {log.log} \
        2> {log.err}
        """


# Generate list of read IDs to keep ("no ribo") from the filtered R2 file
rule ilmn_2a_bwa_rRNA_get_no_rRNA_list:
    input:
        R2_FQ_BWA_FILTERED="{OUTDIR}/{SAMPLE}/short_read/rRNA/bwa/no_rRNA_R2.fq",
    output:
        NORIBO_LIST=temp("{OUTDIR}/{SAMPLE}/short_read/rRNA/bwa/no_rRNA_readID.list"),
    threads: 1
    shell:
        """
        cat {input.R2_FQ_BWA_FILTERED} | awk -f scripts/awk/fq_readHeader.awk - > {output.NORIBO_LIST}
        """


# Filter R1 reads using the generated list of read IDs
rule ilmn_2a_bwa_rRNA_filter_R1:
    input:
        R1_FQ="{OUTDIR}/{SAMPLE}/short_read/tmp/twiceCut_R1.fq.gz",
        NORIBO_LIST="{OUTDIR}/{SAMPLE}/short_read/rRNA/bwa/no_rRNA_readID.list",
    output:
        R1_FQ_NORIBO="{OUTDIR}/{SAMPLE}/short_read/rRNA/bwa/no_rRNA_R1.fq",
    resources:
        mem="64G",
    threads: config["CORES"]
    conda:
        f"{workflow.basedir}/envs/seqkit.yml"
    shell:
        """
        zcat {input.R1_FQ} | seqkit grep -f {input.NORIBO_LIST} -o {output.R1_FQ_NORIBO}
        """


# Filter R1 reads using the generated list of read IDs
rule ilmn_2a_bwa_rRNA_filter_trimmed_R1:
    input:
        R1_FQ="{OUTDIR}/{SAMPLE}/short_read/tmp/twiceCut_{TRIM}_R1.fq.gz",
        rRNA_LIST="{OUTDIR}/{SAMPLE}/short_read/rRNA/bwa/no_rRNA_readID.list",
    output:
        R1_FQ_NORIBO="{OUTDIR}/{SAMPLE}/short_read/rRNA/bwa/no_rRNA_{TRIM}_R1.fq",
    resources:
        mem="64G",
    threads: config["CORES"]
    conda:
        f"{workflow.basedir}/envs/seqkit.yml"
    shell:
        """
        zcat {input.R1_FQ} | seqkit grep -f {input.NORIBO_LIST} -o {output.R1_FQ_NORIBO}
        """


# Compress previous outputs
rule ilmn_2a_bwa_rRNA_compress_unmapped:
    input:
        R1_FQ_BWA_FILTERED="{OUTDIR}/{SAMPLE}/short_read/rRNA/bwa/no_rRNA_{READ}.fq",
    output:
        R1_FQ_BWA_FILTERED="{OUTDIR}/{SAMPLE}/short_read/rRNA/bwa/no_rRNA_{READ}.fq.gz",
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
        FQ="{OUTDIR}/{SAMPLE}/short_read/rRNA/bwa/no_rRNA_{READ}.fq.gz",
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
