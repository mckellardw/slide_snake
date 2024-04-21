# fastqc before trimming
rule ont_nanoQC:
    input:
        MERGED_FQ="{OUTDIR}/{SAMPLE}/tmp/ont/merged.fq.gz",
    output:
        HTML="{OUTDIR}/{SAMPLE}/ont/nanoqc/ont_preAdapterScan/nanoQC.html",
    params:
        MIN_LENGTH=8,
        NANOQC=EXEC["NANOQC"],
    # config['CORES']
    threads: 1
    log:
        log="{OUTDIR}/{SAMPLE}/ont/nanoqc/ont_preAdapterScan/NanoQC.log",
    conda:
        f"{workflow.basedir}/envs/nanoQC.yml"
    shell:
        """
        mkdir -p $(dirname {output.HTML})

        {params.NANOQC} \
            --outdir $(dirname {output.HTML}) \
            --minlen {params.MIN_LENGTH} \
            {input.MERGED_FQ} \
        2>> {log.log}
        """


# nanoqc after cutadapt trimming
rule ont_nanoQC_postCutadapt:
    input:
        FQ="{OUTDIR}/{SAMPLE}/tmp/ont/cut_{READ}.fq.gz",
    output:
        HTML="{OUTDIR}/{SAMPLE}/ont/nanoqc/ont_postCutadapt_{READ}/nanoQC.html",
    params:
        MIN_LENGTH=8,
        NANOQC=EXEC["NANOQC"],
    threads: config["CORES"]
    # min([config['CORES'],8]) # 8 core max
    log:
        log="{OUTDIR}/{SAMPLE}/ont/nanoqc/ont_postCutadapt_{READ}/NanoQC.log",
    conda:
        f"{workflow.basedir}/envs/nanoQC.yml"
    shell:
        """
        mkdir -p $(dirname {output.HTML})

        {params.NANOQC} \
            --outdir $(dirname {output.HTML}) \
            --minlen {params.MIN_LENGTH} \
            {input.FQ} \
        2>> {log.log}
        """
