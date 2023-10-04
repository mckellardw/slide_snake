#########################################################
## kallisto pseudoalignment for RNA velocity analysis
#########################################################
# Source: https://bustools.github.io/BUS_notebooks_R/velocity.html


rule kallisto_align_velocity:
    input:
        R1_FQ = '{OUTDIR}/{sample}/tmp/{sample}_R1_final.fq.gz',
        R2_FQ = '{OUTDIR}/{sample}/tmp/{sample}_R2_final.fq.gz',
        R1_FQ_FILTERED = '{OUTDIR}/{sample}/tmp/{sample}_R1_final_filtered.fq.gz',
        R2_FQ_FILTERED = '{OUTDIR}/{sample}/tmp/{sample}_R2_final_filtered.fq.gz',
        BB = "{OUTDIR}/{sample}/bb/whitelist.txt"
    output:
        BUS = temp('{OUTDIR}/{sample}/kb_velo/output.corrected.bus'),
        TRANSCRIPTS = '{OUTDIR}/{sample}/kb_velo/transcripts.txt',
        ECMAP = temp('{OUTDIR}/{sample}/kb_velo/matrix.ec')
    params:
        MEMLIMIT = config['MEMLIMIT_GB']
    log:
        '{OUTDIR}/{sample}/kb_velo/kallisto_align.log'
    threads:
        config['CORES']
    priority:
        42
    run:
        recipe = RECIPE_DICT[wildcards.sample]
        KB_IDX = IDX_VELO_DICT[wildcards.sample]
        # KB_IDX ="/workdir/dwm269/genomes/mm39_all/kallisto_velo_GRCm39_GENCODEM32_REOT1L-as/mixed.idx" #temp hardcoded
        BB_WHITELIST = f"{input.BB}"
        
        KB_X = RECIPE_SHEET["kb.x"][recipe]

        # Select R2 based on alignment recipe
        if "rRNA" in recipe: # Use trimmed & rRNA-filtered .fq's
            R1 = input.R1_FQ_FILTERED
            R2 = input.R2_FQ_FILTERED
        else: # just trimmed .fq's
            R1 = input.R1_FQ
            R2 = input.R2_FQ

        shell(
            f"""
            bash scripts/kb.sh {OUTDIR}/{wildcards.sample}/kb_velo \
            {KB_IDX} \
            {BB_WHITELIST} \
            {KB_X} \
            {log} \
            {threads} \
            {params.MEMLIMIT} \
            {R1} {R2}
            """
        )

# Split the .bus file for spliced/unspliced outputs
#TODO- split this into two rules for better parallelization in runs
#TODO - fix hardcoded bits...
rule split_bus_velocity:
    input:
        BUS = '{OUTDIR}/{sample}/kb_velo/output.corrected.bus',
        TRANSCRIPTS = '{OUTDIR}/{sample}/kb_velo/transcripts.txt',
        ECMAP = '{OUTDIR}/{sample}/kb_velo/matrix.ec'
    output:
        SPLICED = '{OUTDIR}/{sample}/kb_velo/spliced.bus',
        UNSPLICED = '{OUTDIR}/{sample}/kb_velo/unspliced.bus'
    params:
        MATDIR = directory('{OUTDIR}/{sample}/kb_velo/raw')
    threads:
        1
    run:
        # KB_IDX = IDX_VELO_DICT[wildcards.sample]
        shell(
            f"""
            mkdir -p {params.MATDIR}

            {BUST_EXEC} capture \
            -s -x \
            -o {output.SPLICED} \
            -c /workdir/dwm269/genomes/mm39_all/kallisto_velo_GRCm39_GENCODEM32_REOT1L-as/cDNA.t2c \
            --ecmap {input.ECMAP} \
            --txnames {input.TRANSCRIPTS} \
            {input.BUS}
            """
        )

        shell(
            f"""
            {BUST_EXEC} capture \
            -s -x \
            -o {output.UNSPLICED} \
            -c /workdir/dwm269/genomes/mm39_all/kallisto_velo_GRCm39_GENCODEM32_REOT1L-as/introns.t2c \
            --ecmap {input.ECMAP} \
            --txnames {input.TRANSCRIPTS} \
            {input.BUS}
            """
        )

#
rule bus2mat_velocity:
    input:
        SPLICED = '{OUTDIR}/{sample}/kb_velo/spliced.bus',
        UNSPLICED = '{OUTDIR}/{sample}/kb_velo/unspliced.bus',
        TRANSCRIPTS = '{OUTDIR}/{sample}/kb_velo/transcripts.txt',
        ECMAP = '{OUTDIR}/{sample}/kb_velo/matrix.ec'
    output:
        BCS = '{OUTDIR}/{sample}/kb_velo/raw/output.barcodes.txt',
        GENES = '{OUTDIR}/{sample}/kb_velo/raw/output.genes.txt',
        SPLICED = '{OUTDIR}/{sample}/kb_velo/raw/spliced.mtx',
        UNSPLICED = '{OUTDIR}/{sample}/kb_velo/raw/unspliced.mtx'
    params:
        MATDIR = directory('{OUTDIR}/{sample}/kb_velo/raw')
    threads:
        1
    run:
        shell(
            f"""
            mkdir -p {params.MATDIR}

            {BUST_EXEC} count \
            --output {params.MATDIR}/ \
            --genemap {T2G_VELO_DICT[wildcards.sample]} \
            --ecmap {input.ECMAP} \
            --txnames {input.TRANSCRIPTS} \
            --genecounts \
            --umi-gene \
            --em \
            {input.SPLICED}

            mv {output.SPLICED.replace("spliced","output")} {output.SPLICED}
            """
        )

        shell(
            f"""
            {BUST_EXEC} count \
            --output {params.MATDIR}/ \
            --genemap {T2G_VELO_DICT[wildcards.sample]} \
            --ecmap {input.ECMAP} \
            --txnames {input.TRANSCRIPTS} \
            --genecounts \
            --umi-gene \
            --em \
            {input.UNSPLICED}

            mv {output.UNSPLICED.replace("unspliced","output")} {output.UNSPLICED}
            """
        )

# gzip the count matrix, etc.
rule compress_kb_outs_velocity:
    input:
        BCS = '{OUTDIR}/{sample}/kb_velo/raw/output.barcodes.txt',
        GENES = '{OUTDIR}/{sample}/kb_velo/raw/output.genes.txt',
        SPLICED = '{OUTDIR}/{sample}/kb_velo/raw/spliced.mtx',
        UNSPLICED = '{OUTDIR}/{sample}/kb_velo/raw/unspliced.mtx'
    output:
        BCS = '{OUTDIR}/{sample}/kb_velo/raw/output.barcodes.txt.gz',
        GENES = '{OUTDIR}/{sample}/kb_velo/raw/output.genes.txt.gz',
        SPLICED = '{OUTDIR}/{sample}/kb_velo/raw/spliced.mtx.gz',
        UNSPLICED = '{OUTDIR}/{sample}/kb_velo/raw/unspliced.mtx.gz'
    params:
        MATDIR = directory('{OUTDIR}/{sample}/kb_velo/raw')
    threads:
        config['CORES']        
    run:
        shell(
            f"""
            pigz -p{threads} {input.BCS} {input.GENES} {input.SPLICED} {input.UNSPLICED}
            """
        )