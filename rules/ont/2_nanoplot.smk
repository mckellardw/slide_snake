# fastqc before trimming
rule ont_NanoPlot_preTrim:
    input:
        MERGED_FQ="{OUTDIR}/{SAMPLE}/tmp/ont/merged.fq.gz",
    output:
        DIR=directory("{OUTDIR}/{SAMPLE}/nanoplot/ont_preAdapterScan"),
        SUMMARY="{OUTDIR}/{SAMPLE}/nanoplot/ont_preAdapterScan/summary.txt",
    params:
        adapters=config["FASTQC_ADAPTERS"],
    threads: config["CORES"]
    # min([config['CORES'],8]) # 8 core max
    run:
        shell(
            f"""
            mkdir -p {output.DIR}

            {EXEC['NANOPLOT']} \
                -o {output.DIR} \
                -t {threads} \
                --sumary {output.SUMMARY} \
                --plots dot \
                --dpi 200 \
                --fastq {input.MERGED_FQ}
            """
        )


# fastqc after cutadapt trimming
rule ont_NanoPlot_postCutadapt:
    input:
        FQ="{OUTDIR}/{SAMPLE}/tmp/ont/cut_{READ}.fq.gz",
    output:
        DIR=directory("{OUTDIR}/{SAMPLE}/nanoplot/ont_postCutadapt_{READ}"),
        SUMMARY="{OUTDIR}/{SAMPLE}/nanoplot/ont_postCutadapt_{READ}/summary.txt",
    params:
        adapters=config["FASTQC_ADAPTERS"],
    threads: config["CORES"]
    # min([config['CORES'],8]) # 8 core max
    run:
        shell(
            f"""
            mkdir -p {output.DIR}

            {EXEC['NANOPLOT']} \
                -o {output.DIR} \
                -t {threads} \
                --sumary {output.SUMMARY} \
                --plots dot \
                --dpi 200 \
                --fastq {input.FQ}
            """
        )
