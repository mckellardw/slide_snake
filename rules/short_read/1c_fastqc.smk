#############################################
## fastQC rules
#############################################
# fastqc on R1 & R2 before trimming
rule fastQC_preTrim:
    input:
        MERGED_FQ="{OUTDIR}/{SAMPLE}/tmp/merged_{READ}.fq.gz",
    output:
        fastqcDir=directory("{OUTDIR}/{SAMPLE}/fastqc/preCutadapt_{READ}"),
    params:
        adapters=config["FASTQC_ADAPTERS"],
    resources:
        mem="16G",
        threads=config["CORES"],
    conda:
        f"{workflow.basedir}/envs/fastqc.yml"
    shell:
        """
        mkdir -p {output.fastqcDir}

        fastqc \
            --outdir {output.fastqcDir} \
            --threads {resources.threads} \
            -a {params.adapters} \
            {input.MERGED_FQ}
        """


# fastqc on R1 after linker removal & R2 trimming/filtering
rule fastQC_postTrim:
    input:
        FINAL_FQ="{OUTDIR}/{SAMPLE}/tmp/cut_{READ}.fq.gz",
    output:
        fastqcDir=directory("{OUTDIR}/{SAMPLE}/fastqc/postCutadapt_{READ}"),
    params:
        adapters=config["FASTQC_ADAPTERS"],
    resources:
        mem="16G",
        threads=config["CORES"],
    conda:
        f"{workflow.basedir}/envs/fastqc.yml"
    shell:
        """
        mkdir -p {output.fastqcDir}

        fastqc \
            --outdir {output.fastqcDir} \
            --threads {resources.threads} \
            -a {params.adapters} \
            {input.FINAL_FQ}
        """


rule fastQC_twiceTrim:
    input:
        FINAL_FQ="{OUTDIR}/{SAMPLE}/tmp/twiceCut_{READ}.fq.gz",
    output:
        fastqcDir=directory("{OUTDIR}/{SAMPLE}/fastqc/twiceCutadapt_{READ}"),
    params:
        adapters=config["FASTQC_ADAPTERS"],
    resources:
        mem="16G",
        threads=config["CORES"],
    conda:
        f"{workflow.basedir}/envs/fastqc.yml"
    shell:
        """
        mkdir -p {output.fastqcDir}

        fastqc \
            --outdir {output.fastqcDir} \
            --threads {resources.threads} \
            -a {params.adapters} \
            {input.FINAL_FQ}
        """
