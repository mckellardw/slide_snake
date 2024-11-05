# rRNA removeal w/ RiboDetector
## Source: https://github.com/hzi-bifo/RiboDetector
## Paper: https://doi.org/10.1093/nar/gkac112


# Run ribodetector on preprocessed reads
# TODO- refactor to incorporate internal trimming options into rRNA filtering
rule ilmn_2b_ribodetector:
    input:
        # R1_FQ="{OUTDIR}/{SAMPLE}/short_read/tmp/twiceCut_R1.fq.gz",
        R2_FQ="{OUTDIR}/{SAMPLE}/short_read/tmp/twiceCut_R2.fq.gz",
    output:
        # R1_FQ="{OUTDIR}/{SAMPLE}/short_read/rRNA/ribodetector/noRibo_R1.fq",
        R2_FQ="{OUTDIR}/{SAMPLE}/short_read/rRNA/ribodetector/noRibo_R2.fq",
        # RIBO_R1_FQ="{OUTDIR}/{SAMPLE}/short_read/rRNA/ribodetector/noRibo_R1.fq",
        RIBO_R2_FQ="{OUTDIR}/{SAMPLE}/short_read/rRNA/ribodetector/yesRibo_R2.fq",
    params:
        CHUNK_SIZE=2048,
    log:
        log="{OUTDIR}/{SAMPLE}/short_read/rRNA/ribodetector/ribodetector.log",
    threads: config["CORES"]
    conda:
        f"{workflow.basedir}/envs/ribodetector.yml"
    shell:
        """
        READ_LEN=$(seqkit stats {input.R2_FQ} | tail -n 1 | awk '{{print $7}}')

        ribodetector_cpu \
            --threads {threads} \
            --len  $(printf "%.0f" $READ_LEN) \
            --input {input.R2_FQ} \
            --ensure rrna \
            --chunk_size {params.CHUNK_SIZE} \
            --rrna {output.RIBO_R2_FQ} \
            --output {output.R2_FQ} \
            --log {log.log}
        """


# Get list of read IDs to keep ("no ribo") from the output R2 file
rule ilmn_2b_ribodetector_get_noRibo_list:
    input:
        R2_FQ_NORIBO="{OUTDIR}/{SAMPLE}/short_read/rRNA/ribodetector/noRibo_R2.fq",
    output:
        NORIBO_LIST=temp("{OUTDIR}/{SAMPLE}/short_read/rRNA/ribodetector/noRibo_readID.list"),
    threads: 1
    shell:
        """
        cat {input.R2_FQ_NORIBO} \
        | awk -f scripts/awk/fq_readIDs.awk - \
        > {output.NORIBO_LIST}
        """


# Temporarily decompress R1 fastq...
rule ilmn_2b_ribodetector_gunzip_R1:
    input:
        R1_FQ="{OUTDIR}/{SAMPLE}/short_read/tmp/{TMP}_R1.fq.gz",
    output:
        R1_FQ=temp("{OUTDIR}/{SAMPLE}/short_read/tmp/{TMP}_R1.fq"),
    threads: 1
    shell:
        """
        zcat {input.R1_FQ} > {output.R1_FQ} 
        """


# Keep only reads which do NOT contain rRNA sequences
rule ilmn_2b_ribodetector_filter_R1:
    input:
        R1_FQ="{OUTDIR}/{SAMPLE}/short_read/tmp/twiceCut_R1.fq",
        NORIBO_LIST="{OUTDIR}/{SAMPLE}/short_read/rRNA/ribodetector/noRibo_readID.list",
    output:
        R1_FQ_NORIBO="{OUTDIR}/{SAMPLE}/short_read/rRNA/ribodetector/noRibo_R1.fq",
    resources:
        mem="32G",
    threads: 1
    shell:
        """
        seqtk subseq {input.R1_FQ} {input.NORIBO_LIST} \
        > {output.R1_FQ_NORIBO}
        """


# Internally trimmed reads
rule ilmn_2b_ribodetector_filter_R1_internalTrim:
    input:
        R1_FQ="{OUTDIR}/{SAMPLE}/short_read/tmp/twiceCut_internalTrim_R1.fq",
        NORIBO_LIST="{OUTDIR}/{SAMPLE}/short_read/rRNA/ribodetector/noRibo_readID.list",
    output:
        R1_FQ_NORIBO="{OUTDIR}/{SAMPLE}/short_read/rRNA/ribodetector/noRibo_internalTrim_R1.fq",
    resources:
        mem="32G",
    threads: 1
    shell:
        """
        seqtk subseq {input.R1_FQ} {input.NORIBO_LIST} \
        > {output.R1_FQ_NORIBO}
        """


# Hard-trimmed reads
rule ilmn_2b_ribodetector_filter_R1_hardTrim:
    input:
        R1_FQ="{OUTDIR}/{SAMPLE}/short_read/tmp/twiceCut_hardTrim_R1.fq",
        NORIBO_LIST="{OUTDIR}/{SAMPLE}/short_read/rRNA/ribodetector/noRibo_readID.list",
    output:
        R1_FQ_NORIBO="{OUTDIR}/{SAMPLE}/short_read/rRNA/ribodetector/noRibo_hardTrim_R1.fq",
    resources:
        mem="32G",
    threads: 1
    shell:
        """
        seqtk subseq {input.R1_FQ} {input.NORIBO_LIST} \
        > {output.R1_FQ_NORIBO}
        """


# Compress all the R1 files
rule ilmn_2b_ribodetector_compress_fqs:
    input:
        FQ="{OUTDIR}/{SAMPLE}/short_read/rRNA/ribodetector/noRibo_{READ}.fq",
    output:
        FQ="{OUTDIR}/{SAMPLE}/short_read/rRNA/ribodetector/noRibo_{READ}.fq.gz",
    resources:
        mem="8G",
    threads: config["CORES"]
    shell:
        """
        pigz -p{threads} {input.FQ}
        """
