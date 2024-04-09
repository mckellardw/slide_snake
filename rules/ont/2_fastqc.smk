# fastqc before trimming
rule ont_fastQC_preTrim:
    input:
        MERGED_FQ="{OUTDIR}/{SAMPLE}/tmp/ont/merged.fq.gz",
    output:
        DIR=directory("{OUTDIR}/{SAMPLE}/fastqc/ont_preAdapterScan"),
    params:
        adapters=config["FASTQC_ADAPTERS"],
    threads: config["CORES"]
    # min([config['CORES'],8]) # 8 core max
    run:
        shell(
            f"""
            mkdir -p {output.DIR}

            {EXEC['FASTQC']} \
                --outdir {output.DIR} \
                --threads {threads} \
                -a {params.adapters} \
                {input.MERGED_FQ}
            """
        )


# fastqc after cutadapt trimming
rule ont_fastQC_postCutadapt:
    input:
        MERGED_FQ="{OUTDIR}/{SAMPLE}/tmp/ont/cut_{READ}.fq.gz",
    output:
        DIR=directory("{OUTDIR}/{SAMPLE}/fastqc/ont_postCutadapt_{READ}"),
    params:
        adapters=config["FASTQC_ADAPTERS"],
    threads: config["CORES"]
    # min([config['CORES'],8]) # 8 core max
    run:
        shell(
            f"""
            mkdir -p {output.DIR}

            {EXEC['FASTQC']} \
                --outdir {output.DIR} \
                --threads {threads} \
                -a {params.adapters} \
                {input.MERGED_FQ}
            """
        )
