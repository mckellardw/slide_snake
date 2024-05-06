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
    threads: config["CORES"]
    # min([config['CORES'],8]) # 8 core max
    run:
        shell(
            f"""
            mkdir -p {output.fastqcDir}

            {EXEC['FASTQC']} \
                --outdir {output.fastqcDir} \
                --threads {threads} \
                -a {params.adapters} \
                {input.MERGED_FQ}
            """
        )


# fastqc on R1 after linker removal & R2 trimming/filtering
rule fastQC_postTrim:
    input:
        FINAL_FQ="{OUTDIR}/{SAMPLE}/tmp/cut_{READ}.fq.gz",
    output:
        fastqcDir=directory("{OUTDIR}/{SAMPLE}/fastqc/postCutadapt_{READ}"),
        # fastqcReport = ''
    threads: config["CORES"]
    # min([config['CORES'],8]) # 8 core max
    params:
        adapters=config["FASTQC_ADAPTERS"],
    run:
        shell(
            f"""
            mkdir -p {output.fastqcDir}

            {EXEC['FASTQC']} \
                --outdir {output.fastqcDir} \
                --threads {threads} \
                -a {params.adapters} \
                {input.FINAL_FQ}
            """
        )


rule fastQC_twiceTrim:
    input:
        FINAL_FQ="{OUTDIR}/{SAMPLE}/tmp/twiceCut_{READ}.fq.gz",
    output:
        fastqcDir=directory("{OUTDIR}/{SAMPLE}/fastqc/twiceCutadapt_{READ}"),
        # fastqcReport = ''
    threads: config["CORES"]
    # min([config['CORES'],8]) # 8 core max
    params:
        adapters=config["FASTQC_ADAPTERS"],
    run:
        shell(
            f"""
            mkdir -p {output.fastqcDir}

            {EXEC['FASTQC']} \
                --outdir {output.fastqcDir} \
                --threads {threads} \
                -a {params.adapters} \
                {input.FINAL_FQ}
            """
        )