# fastqc before trimming
rule ont_fastQC_preTrim:
    input:
        FQ="{OUTDIR}/{SAMPLE}/tmp/ont/merged.fq.gz",
    output:
        DIR=directory("{OUTDIR}/{SAMPLE}/fastqc/ont_preAdapterScan"),
    params:
        adapters=config["FASTQC_ADAPTERS"],
    threads: config["CORES"]
    # min([config['CORES'],8]) # 8 core max
    conda:
        f"{workflow.basedir}/envs/fastqc.yml"
    shell:
        """
        mkdir -p {output.DIR}

        fastqc \
            --outdir {output.DIR} \
            --threads {threads} \
            -a {params.adapters} \
            {input.FQ}
        """


# fastqc after cutadapt trimming
rule ont_fastQC_preCutadapt:
    input:
        FQ="{OUTDIR}/{SAMPLE}/tmp/ont/adapter_scan_readids/merged_adapter_{READ}.fq.gz",
    output:
        DIR=directory("{OUTDIR}/{SAMPLE}/fastqc/ont_preCutadapt_{READ}"),
    params:
        adapters=config["FASTQC_ADAPTERS"],
    threads: config["CORES"]
    # min([config['CORES'],8]) # 8 core max
    conda:
        f"{workflow.basedir}/envs/fastqc.yml"
    shell:
        """
        mkdir -p {output.DIR}

        fastqc \
            --outdir {output.DIR} \
            --threads {threads} \
            -a {params.adapters} \
            {input.FQ}
        """


# fastqc after cutadapt trimming
rule ont_fastQC_postCutadapt:
    input:
        FQ="{OUTDIR}/{SAMPLE}/tmp/ont/cut_{READ}.fq.gz",
    output:
        DIR=directory("{OUTDIR}/{SAMPLE}/fastqc/ont_postCutadapt_{READ}"),
    params:
        adapters=config["FASTQC_ADAPTERS"],
    threads: config["CORES"]
    # min([config['CORES'],8]) # 8 core max
    conda:
        f"{workflow.basedir}/envs/fastqc.yml"
    shell:
        """
            mkdir -p {output.DIR}

            fastqc \
                --outdir {output.DIR} \
                --threads {threads} \
                -a {params.adapters} \
                {input.FQ}
            """
        