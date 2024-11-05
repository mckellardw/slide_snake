# fastqc before trimming
rule ont_fastQC_preTrim:
    input:
        FQ="{OUTDIR}/{SAMPLE}/ont/tmp/merged.fq.gz",
    output:
        DIR=directory("{OUTDIR}/{SAMPLE}/ont/fastqc/ont_preAdapterScan"),
    params:
        adapters=config["FASTQC_ADAPTERS"],
    threads: config["CORES"]
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
        FQ="{OUTDIR}/{SAMPLE}/ont/tmp/adapter_scan_readids/merged_adapter_{READ}.fq.gz",
    output:
        DIR=directory("{OUTDIR}/{SAMPLE}/ont/fastqc/ont_preCutadapt_{READ}"),
    params:
        adapters=config["FASTQC_ADAPTERS"],
    threads: config["CORES"]
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
        FQ="{OUTDIR}/{SAMPLE}/ont/tmp/cut_{READ}.fq.gz",
    output:
        DIR=directory("{OUTDIR}/{SAMPLE}/ont/fastqc/ont_postCutadapt_{READ}"),
    params:
        adapters=config["FASTQC_ADAPTERS"],
    threads: config["CORES"]
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
