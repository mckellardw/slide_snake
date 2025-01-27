# rRNA removeal w/ RiboDetector
## Source: https://github.com/hzi-bifo/RiboDetector
## Paper: https://doi.org/10.1093/nar/gkac112


# Run ribodetector on preprocessed reads
# TODO- refactor to incorporate internal trimming options into rRNA filtering
rule ilmn_2b_ribodetector:
    input:
        R2_FQ="{OUTDIR}/{SAMPLE}/short_read/tmp/twiceCut_R2.fq.gz",
    output:
        R2_FQ="{OUTDIR}/{SAMPLE}/short_read/rRNA/ribodetector/noRibo_R2.fq",
        RIBO_R2_FQ="{OUTDIR}/{SAMPLE}/short_read/rRNA/ribodetector/yesRibo_R2.fq",
    params:
        CHUNK_SIZE=1024,
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
        NORIBO_LIST=temp(
            "{OUTDIR}/{SAMPLE}/short_read/rRNA/ribodetector/noRibo_readID.list"
        ),
    threads: 1
    shell:
        """
        cat {input.R2_FQ_NORIBO} \
        | awk -f scripts/awk/fq_readIDs.awk - \
        > {output.NORIBO_LIST}
        """


# Combine filter R1 reads with wildcard
rule ilmn_2b_ribodetector_filter_R1:
    input:
        R1_FQ="{OUTDIR}/{SAMPLE}/short_read/tmp/twiceCut_R1.fq.gz",
        NORIBO_LIST="{OUTDIR}/{SAMPLE}/short_read/rRNA/ribodetector/noRibo_readID.list",
    output:
        R1_FQ_NORIBO="{OUTDIR}/{SAMPLE}/short_read/rRNA/ribodetector/noRibo_R1.fq",
    resources:
        mem="32G",
    threads: 1
    conda:
        f"{workflow.basedir}/envs/seqkit.yml"
    shell:
        """
        zcat {input.R1_FQ} | seqkit grep -f {input.NORIBO_LIST} -o {output.R1_FQ_NORIBO}
        """


rule ilmn_2b_ribodetector_filter_trimmed_R1:
    input:
        R1_FQ="{OUTDIR}/{SAMPLE}/short_read/tmp/twiceCut_{TRIM}_R1.fq.gz",
        NORIBO_LIST="{OUTDIR}/{SAMPLE}/short_read/rRNA/ribodetector/noRibo_readID.list",
    output:
        R1_FQ_NORIBO="{OUTDIR}/{SAMPLE}/short_read/rRNA/ribodetector/noRibo_{TRIM}_R1.fq",
    resources:
        mem="32G",
    threads: 1
    conda:
        f"{workflow.basedir}/envs/seqkit.yml"
    shell:
        """
        zcat {input.R1_FQ} | seqkit grep -f {input.NORIBO_LIST} -o {output.R1_FQ_NORIBO}
        """


# Compress all the ribodetector fastqs
rule ilmn_2b_ribodetector_compress_fqs:
    input:
        FQ="{OUTDIR}/{SAMPLE}/short_read/rRNA/ribodetector/{RIBO}Ribo_{READ}.fq",
    output:
        FQ="{OUTDIR}/{SAMPLE}/short_read/rRNA/ribodetector/{RIBO}Ribo_{READ}.fq.gz",
    resources:
        mem="8G",
    threads: config["CORES"]
    shell:
        """
        pigz -p{threads} {input.FQ}
        """
