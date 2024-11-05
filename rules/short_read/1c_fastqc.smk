#############################################
## fastQC rules
#############################################
# fastqc on R1 & R2 before trimming
rule ilmn_1c_fastQC_preTrim:
    input:
        MERGED_FQ="{OUTDIR}/{SAMPLE}/short_read/tmp/merged_{READ}.fq.gz",
    output:
        FASTQC_DIR=directory("{OUTDIR}/{SAMPLE}/short_read/fastqc/preCutadapt_{READ}"),
    params:
        adapters=config["FASTQC_ADAPTERS"],
    resources:
        mem="16G",
    threads: config["CORES"]
    conda:
        f"{workflow.basedir}/envs/fastqc.yml"
    shell:
        """
        mkdir -p {output.FASTQC_DIR}

        fastqc \
            --outdir {output.FASTQC_DIR} \
            --threads {threads} \
            -a {params.adapters} \
            {input.MERGED_FQ}
        """


# fastqc on R1 after linker removal & R2 trimming/filtering
rule ilmn_1c_fastQC_postTrim:
    input:
        FINAL_FQ="{OUTDIR}/{SAMPLE}/short_read/tmp/cut_{READ}.fq.gz",
    output:
        FASTQC_DIR=directory("{OUTDIR}/{SAMPLE}/short_read/fastqc/postCutadapt_{READ}"),
    params:
        adapters=config["FASTQC_ADAPTERS"],
    resources:
        mem="16G",
    threads: config["CORES"]
    conda:
        f"{workflow.basedir}/envs/fastqc.yml"
    shell:
        """
        mkdir -p {output.FASTQC_DIR}

        fastqc \
            --outdir {output.FASTQC_DIR} \
            --threads {threads} \
            -a {params.adapters} \
            {input.FINAL_FQ}
        """


rule ilmn_1c_fastQC_twiceTrim:
    input:
        FINAL_FQ="{OUTDIR}/{SAMPLE}/short_read/tmp/twiceCut_{READ}.fq.gz",
    output:
        FASTQC_DIR=directory("{OUTDIR}/{SAMPLE}/short_read/fastqc/twiceCutadapt_{READ}"),
    params:
        adapters=config["FASTQC_ADAPTERS"],
    resources:
        mem="16G",
    threads: config["CORES"]
    conda:
        f"{workflow.basedir}/envs/fastqc.yml"
    shell:
        """
        mkdir -p {output.FASTQC_DIR}

        fastqc \
            --outdir {output.FASTQC_DIR} \
            --threads {threads} \
            -a {params.adapters} \
            {input.FINAL_FQ}
        """
