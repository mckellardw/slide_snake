#############################################
## Read trimming & QC rules
#############################################
# fastqc before trimming
rule ont_fastQC_preTrim:
    input:
        MERGED_FQ = '{OUTDIR}/{SAMPLE}/ont/merged.fq.gz'
    output:
        DIR = directory('{OUTDIR}/{SAMPLE}/fastqc/ont_preAdapterScan')
    params:
        adapters = config['FASTQC_ADAPTERS']
    threads:
        config['CORES']
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

# https://github.com/yfukasawa/LongQC
rule ont_longQC_preTrim:
    input:
        MERGED_FQ = '{OUTDIR}/{SAMPLE}/ont/merged.fq.gz'
    output:
        DIR = directory('{OUTDIR}/{SAMPLE}/longqc/ont_preAdapterScan')
    params:
        adapters = config['FASTQC_ADAPTERS']
    threads:
        config['CORES']
        # min([config['CORES'],8]) # 8 core max
    run:
        shell(
            f"""
            mkdir -p {output.DIR}

            {EXEC['LONGQC']} \
                -o {output.DIR} \
                -p {threads} \
                -x ont-ligation \
                {input.MERGED_FQ}
            """
        )