#############################################
## kallisto pseudoalignment
#############################################

rule kallisto_align:
    input:
        R1_FQ = '{OUTDIR}/{SAMPLE}/tmp/final_R1.fq.gz',
        R2_FQ = '{OUTDIR}/{SAMPLE}/tmp/final_R2.fq.gz',
        R1_FQ_FILTERED = '{OUTDIR}/{SAMPLE}/tmp/final_filtered_R1.fq.gz',
        R2_FQ_FILTERED = '{OUTDIR}/{SAMPLE}/tmp/final_filtered_R2.fq.gz',
        BB = "{OUTDIR}/{SAMPLE}/bb/whitelist.txt"
    output:
        BUS = temp('{OUTDIR}/{SAMPLE}/kb/{RECIPE}/output.bus'),
        BUS_CORRECTED = temp('{OUTDIR}/{SAMPLE}/kb/{RECIPE}/output.corrected.bus'),
        TRANSCRIPTS = '{OUTDIR}/{SAMPLE}/kb/{RECIPE}/transcripts.txt',
        ECMAP = temp('{OUTDIR}/{SAMPLE}/kb/{RECIPE}/matrix.ec')
    params:
        MEMLIMIT = config['MEMLIMIT_GB']
    log:
        '{OUTDIR}/{SAMPLE}/kb/{RECIPE}/kallisto_align.log'
    threads:
        config['CORES']
    resources:
        mem_mb = config['MEMLIMIT_MB']
    priority:
        42
    run:
        tmp_recipe = RECIPE_DICT[wildcards.SAMPLE]
        KB_IDX = IDX_DICT[wildcards.SAMPLE]
        
        KB_X = RECIPE_SHEET["kb.x"][tmp_recipe]

        # Select R2 based on alignment recipe
        if "rRNA" in tmp_recipe: # Use trimmed & rRNA-filtered .fq's
            R1 = input.R1_FQ_FILTERED
            R2 = input.R2_FQ_FILTERED
        else: # just trimmed .fq's
            R1 = input.R1_FQ
            R2 = input.R2_FQ

        shell(
            f"""
            bash scripts/bash/kb.sh \
            --outdir {OUTDIR}/{wildcards.SAMPLE}/kb \
            --kb_idx {KB_IDX} \
            --whitelist {input.BB} \
            --chemistry {KB_X} \
            --log {log} \
            --threads {threads} \
            --memlimit {params.MEMLIMIT} \
            --r1fq {R1} \
            --r2fq {R2}
            """
        )

rule bus2mat:
    input:
        BUS = '{OUTDIR}/{SAMPLE}/kb/{RECIPE}/output.corrected.bus',
        TRANSCRIPTS = '{OUTDIR}/{SAMPLE}/kb/{RECIPE}/transcripts.txt',
        ECMAP = '{OUTDIR}/{SAMPLE}/kb/{RECIPE}/matrix.ec'
    output:
        BCS = '{OUTDIR}/{SAMPLE}/kb/{RECIPE}/raw/output.barcodes.txt',
        GENES = '{OUTDIR}/{SAMPLE}/kb/{RECIPE}/raw/output.genes.txt',
        MAT = '{OUTDIR}/{SAMPLE}/kb/{RECIPE}/raw/output.mtx'
    params:
        MATDIR = directory('{OUTDIR}/{SAMPLE}/kb/{RECIPE}/raw')
    threads:
        1
    run:
        KB_T2G = T2G_DICT[wildcards.SAMPLE]

        shell(
            f"""
            mkdir -p {params.MATDIR}

            {EXEC['BUSTOOLS']} count \
                --output {params.MATDIR} \
                --genemap {KB_T2G} \
                --ecmap {input.ECMAP} \
                --txnames {input.TRANSCRIPTS} \
                --genecounts \
                --umi-gene \
                --em \
                {input.BUS}
            """
        )

# gzip the count matrix, etc.
rule compress_kb_outs:
    input:
        BCS = '{OUTDIR}/{SAMPLE}/kb/{RECIPE}/raw/output.barcodes.txt',
        GENES = '{OUTDIR}/{SAMPLE}/kb/{RECIPE}/raw/output.genes.txt',
        MAT = '{OUTDIR}/{SAMPLE}/kb/{RECIPE}/raw/output.mtx'
    output:
        BCS = '{OUTDIR}/{SAMPLE}/kb/{RECIPE}/raw/output.barcodes.txt.gz',
        GENES = '{OUTDIR}/{SAMPLE}/kb/{RECIPE}/raw/output.genes.txt.gz',
        MAT = '{OUTDIR}/{SAMPLE}/kb/{RECIPE}/raw/output.mtx.gz'
    params:
        MATDIR = directory('{OUTDIR}/{SAMPLE}/kb/{RECIPE}/raw')
    threads:
        config['CORES']        
    run:
        shell(
            f"""
            {EXEC['PIGZ']} -p{threads} {input.BCS} {input.GENES} {input.MAT}
            """
        )