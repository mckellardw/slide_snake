# Extract rRNA sequences from cDNA fasta file
rule ilmn_2a_extract_rRNA_fasta:
    input:
        CDNA_FA=lambda w: SAMPLE_SHEET["cdna_fa"][w.SAMPLE],
    output:
        FASTA_rRNA="{OUTDIR}/{SAMPLE}/short_read/rRNA/bwa/ref/rRNA_sequences.fa.gz",
    params:
        RRNA_KEYWORDS=config.get(
            "RRNA_KEYWORDS",
            "rRNA,Mt_rRNA,ribosomal_RNA,5S_rRNA,5.8S_rRNA,18S_rRNA,28S_rRNA,12S_rRNA,16S_rRNA",
        ),
    log:
        log="{OUTDIR}/{SAMPLE}/short_read/rRNA/bwa/extract_rRNA.log",
        err="{OUTDIR}/{SAMPLE}/short_read/rRNA/bwa/extract_rRNA.err",
    resources:
        mem="16G",
        time="0:30:00",
    threads: 1
    shell:
        """
        bash scripts/bash/fa_extract_rRNA.sh \
            {input.CDNA_FA} \
            {output.FASTA_rRNA} \
            "{params.RRNA_KEYWORDS}" \
        1> {log.log} \
        2> {log.err}
        """


# Build GTF file from rRNA fasta sequences
rule ilmn_2a_build_rRNA_gtf:
    input:
        CDNA_FA=lambda w: SAMPLE_SHEET["cdna_fa"][w.SAMPLE],
    output:
        GTF_rRNA="{OUTDIR}/{SAMPLE}/short_read/rRNA/bwa/ref/rRNA_annotations.gtf.gz",
    params:
        RRNA_KEYWORDS=config.get(
            "RRNA_KEYWORDS",
            "rRNA,Mt_rRNA,ribosomal_RNA,5S_rRNA,5.8S_rRNA,18S_rRNA,28S_rRNA,12S_rRNA,16S_rRNA",
        ),
    log:
        log="{OUTDIR}/{SAMPLE}/short_read/rRNA/bwa/build_rRNA_gtf.log",
        err="{OUTDIR}/{SAMPLE}/short_read/rRNA/bwa/build_rRNA_gtf.err",
    resources:
        mem="16G",
        time="0:30:00",
    threads: 1
    shell:
        """
        bash scripts/bash/fa_build_rRNA_gtf.sh \
            {input.CDNA_FA} \
            {output.GTF_rRNA} \
            "{params.RRNA_KEYWORDS}" \
        1> {log.log} \
        2> {log.err}
        """


# Build BWA index for rRNA sequences
rule ilmn_2a_build_rRNA_bwa_index:
    input:
        FASTA_rRNA="{OUTDIR}/{SAMPLE}/short_read/rRNA/bwa/ref/rRNA_sequences.fa.gz",
    output:
        BWA_INDEX_FILES=multiext(
            "{OUTDIR}/{SAMPLE}/short_read/rRNA/bwa/ref/rRNA_sequences.fa.gz",
            ".0123",
            ".amb",
            ".ann",
            ".bwt.2bit.64",
            ".pac",
        ),
    log:
        log="{OUTDIR}/{SAMPLE}/short_read/rRNA/bwa/build_bwa_index.log",
        err="{OUTDIR}/{SAMPLE}/short_read/rRNA/bwa/build_bwa_index.err",
    resources:
        mem="32G",
    threads: config["CORES"]
    conda:
        f"{workflow.basedir}/envs/bwa.yml"
    shell:
        """
        bwa-mem2 index {input.FASTA_rRNA} \
        1> {log.log} \
        2> {log.err}
        """


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
        BWA_INDEX_FILES=multiext(
            "{OUTDIR}/{SAMPLE}/short_read/rRNA/bwa/ref/rRNA_sequences.fa.gz",
            ".0123",
            ".amb",
            ".ann",
            ".bwt.2bit.64",
            ".pac",
        ),
        BWA_REF="{OUTDIR}/{SAMPLE}/short_read/rRNA/bwa/ref/rRNA_sequences.fa.gz",
    output:
        BAM1=temp("{OUTDIR}/{SAMPLE}/short_read/rRNA/bwa/aligned.sam"),
        BAM2="{OUTDIR}/{SAMPLE}/short_read/rRNA/bwa/aligned_sorted.bam",
        R2_FQ_BWA_FILTERED="{OUTDIR}/{SAMPLE}/short_read/rRNA/bwa/no_rRNA_R2.fq",
    params:
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
            -f {input.BWA_REF} \
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
        NORIBO_LIST="{OUTDIR}/{SAMPLE}/short_read/rRNA/bwa/no_rRNA_readID.list",
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
