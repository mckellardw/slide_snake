#############################################
## kallisto pseudoalignment
#############################################

rule kallisto_align:
    input:
        FINAL_R1_FQ = '{OUTDIR}/{sample}/tmp/{sample}_R1_final.fq.gz',
        FINAL_R2_FQ = '{OUTDIR}/{sample}/tmp/{sample}_R2_final.fq.gz',
        BB = "{OUTDIR}/{sample}/bb/whitelist.txt"
    output:
        BUSTEXT = temp('{OUTDIR}/{sample}/kb/output.corrected.bus'),
        TRANSCRIPTS = '{OUTDIR}/{sample}/kb/transcripts.txt',
        ECMAP = temp('{OUTDIR}/{sample}/kb/matrix.ec')
    params:
        MEMLIMIT = config['MEMLIMIT']
    log:
        '{OUTDIR}/{sample}/kb/kallisto_align.log'
    threads:
        config['CORES']
    priority:
        42
    run:
        tmp_chemistry = CHEM_DICT[wildcards.sample]
        KB_IDX = IDX_DICT[wildcards.sample]
        BB_WHITELIST = f"{input.BB}"

        KB_X = CHEMISTRY_SHEET["kb.x"][tmp_chemistry]

        shell(
            f"""
            bash scripts/kb.sh {OUTDIR}/{wildcards.sample}/kb \
            {KB_IDX} \
            {BB_WHITELIST} \
            {KB_X} \
            {log} \
            {threads} \
            {params.MEMLIMIT} \
            {input.FINAL_R1_FQ} {input.FINAL_R2_FQ}
            """
        )

rule bus2mat:
    input:
        BUS = '{OUTDIR}/{sample}/kb/output.corrected.bus',
        TRANSCRIPTS = '{OUTDIR}/{sample}/kb/transcripts.txt',
        ECMAP = '{OUTDIR}/{sample}/kb/matrix.ec'
    output:
        MAT = '{OUTDIR}/{sample}/kb/counts_unfiltered/output.mtx'
        # EC = '{OUTDIR}/{sample}/kb/counts_unfiltered/output.ec.txt'
    params:
        MATDIR = directory('{OUTDIR}/{sample}/kb/counts_unfiltered'),
        BUST_EXEC = config['BUST_EXEC']
    threads:
        1
    run:
        KB_T2G = T2G_DICT[wildcards.sample]

        shell(
            f"""
            mkdir -p {params.MATDIR}

            {params.BUST_EXEC} count \
            --output {params.MATDIR}/ \
            --genemap {KB_T2G} \
            --ecmap {input.ECMAP} \
            --txnames {input.TRANSCRIPTS} \
            --genecounts \
            --umi-gene \
            --em \
            {input.BUS}
            """
        )
