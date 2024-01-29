#############################################
## kallisto pseudoalignment
#############################################

rule kallisto_align:
    input:
        R1_FQ = '{OUTDIR}/{sample}/tmp/final_R1.fq.gz',
        R2_FQ = '{OUTDIR}/{sample}/tmp/final_R2.fq.gz',
        R1_FQ_FILTERED = '{OUTDIR}/{sample}/tmp/final_filtered_R1.fq.gz',
        R2_FQ_FILTERED = '{OUTDIR}/{sample}/tmp/final_filtered_R2.fq.gz',
        BB = "{OUTDIR}/{sample}/bb/whitelist.txt"
    output:
        BUS = temp('{OUTDIR}/{sample}/kb/output.bus'),
        BUS_CORRECTED = temp('{OUTDIR}/{sample}/kb/output.corrected.bus'),
        TRANSCRIPTS = '{OUTDIR}/{sample}/kb/transcripts.txt',
        ECMAP = temp('{OUTDIR}/{sample}/kb/matrix.ec')
    params:
        MEMLIMIT = config['MEMLIMIT_GB']
    log:
        '{OUTDIR}/{sample}/kb/kallisto_align.log'
    threads:
        config['CORES']
    resources:
        mem_mb = config['MEMLIMIT_MB']
    priority:
        42
    run:
        tmp_recipe = RECIPE_DICT[wildcards.sample]
        KB_IDX = IDX_DICT[wildcards.sample]
        
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
                --outdir {OUTDIR}/{wildcards.sample}/kb \
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
        BUS = '{OUTDIR}/{sample}/kb/output.corrected.bus',
        TRANSCRIPTS = '{OUTDIR}/{sample}/kb/transcripts.txt',
        ECMAP = '{OUTDIR}/{sample}/kb/matrix.ec'
    output:
        BCS = '{OUTDIR}/{sample}/kb/raw/output.barcodes.txt',
        GENES = '{OUTDIR}/{sample}/kb/raw/output.genes.txt',
        MAT = '{OUTDIR}/{sample}/kb/raw/output.mtx'
        # EC = '{OUTDIR}/{sample}/kb/raw/output.ec.txt'
    params:
        MATDIR = directory('{OUTDIR}/{sample}/kb/raw')
    threads:
        1
    run:
        KB_T2G = T2G_DICT[wildcards.sample]

        shell(
            f"""
            mkdir -p {params.MATDIR}

            {EXEC['BUSTOOLS']} count \
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

# gzip the count matrix, etc.
rule compress_kb_outs:
    input:
        BCS = '{OUTDIR}/{sample}/kb/raw/output.barcodes.txt',
        GENES = '{OUTDIR}/{sample}/kb/raw/output.genes.txt',
        MAT = '{OUTDIR}/{sample}/kb/raw/output.mtx'
    output:
        BCS = '{OUTDIR}/{sample}/kb/raw/output.barcodes.txt.gz',
        GENES = '{OUTDIR}/{sample}/kb/raw/output.genes.txt.gz',
        MAT = '{OUTDIR}/{sample}/kb/raw/output.mtx.gz'
    params:
        MATDIR = directory('{OUTDIR}/{sample}/kb/raw')
    threads:
        config['CORES']        
    run:
        shell(
            f"""
            {EXEC['PIGZ']} -p{threads} {input.BCS} {input.GENES} {input.MAT}
            """
        )