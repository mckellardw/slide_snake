#########################################################
## kallisto pseudoalignment for RNA velocity analysis
#########################################################
# Source: https://bustools.github.io/BUS_notebooks_R/velocity.html
#BUStools manual: https://bustools.github.io/manual

rule kallisto_align_velocity:
    input:
        R1_FQ = '{OUTDIR}/{sample}/tmp/final_R1.fq.gz',
        R2_FQ = '{OUTDIR}/{sample}/tmp/final_R2.fq.gz',
        R1_FQ_FILTERED = '{OUTDIR}/{sample}/tmp/final_filtered_R1.fq.gz',
        R2_FQ_FILTERED = '{OUTDIR}/{sample}/tmp/final_filtered_R2.fq.gz',
        BB = "{OUTDIR}/{sample}/bb/whitelist.txt"
    output:
        BUS = temp('{OUTDIR}/{sample}/kb_velo/output.bus'),
        BUS_CORRECTED = temp('{OUTDIR}/{sample}/kb_velo/output.corrected.bus'),
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
#TODO - fix hardcoded bits...
rule split_bus_velocity_spliced:
    input:
        BUS = '{OUTDIR}/{sample}/kb_velo/output.corrected.bus',
        TRANSCRIPTS = '{OUTDIR}/{sample}/kb_velo/transcripts.txt',
        ECMAP = '{OUTDIR}/{sample}/kb_velo/matrix.ec'
    output:
        SPLICED = '{OUTDIR}/{sample}/kb_velo/spliced.bus'
    params:
        MATDIR = directory('{OUTDIR}/{sample}/kb_velo/raw')
    log:
        '{OUTDIR}/{sample}/kb_velo/split_bus_velocity_spliced.log'
    threads:
        1
    run:
        KB_IDX = IDX_VELO_DICT[wildcards.sample]
        shell(
            f"""
            mkdir -p {params.MATDIR}

            {EXEC['BUSTOOLS']} capture \
                --transcripts \
                --output {output.SPLICED} \
                --capture {KB_IDX}/cDNA.t2c \
                --ecmap {input.ECMAP} \
                --txnames {input.TRANSCRIPTS} \
                {input.BUS}
            """
        )

# Split the .bus file for spliced/unspliced outputs
#TODO - fix hardcoded bits...
rule split_bus_velocity_unspliced:
    input:
        BUS = '{OUTDIR}/{sample}/kb_velo/output.corrected.bus',
        TRANSCRIPTS = '{OUTDIR}/{sample}/kb_velo/transcripts.txt',
        ECMAP = '{OUTDIR}/{sample}/kb_velo/matrix.ec'
    output:
        UNSPLICED = '{OUTDIR}/{sample}/kb_velo/unspliced.bus'
    params:
        MATDIR = directory('{OUTDIR}/{sample}/kb_velo/raw')
    threads:
        1
    run:
        KB_IDX = IDX_VELO_DICT[wildcards.sample]
        shell(
            f"""
            mkdir -p {params.MATDIR}

            {EXEC['BUSTOOLS']} capture \
                --transcripts \
                --output {output.UNSPLICED} \
                --capture {KB_IDX}/introns.t2c \
                --ecmap {input.ECMAP} \
                --txnames {input.TRANSCRIPTS} \
                {input.BUS}
            """
        )

# build spliced matrix
rule bus2mat_velocity_spliced:
    input:
        SPLICED = '{OUTDIR}/{sample}/kb_velo/spliced.bus',
        TRANSCRIPTS = '{OUTDIR}/{sample}/kb_velo/transcripts.txt',
        ECMAP = '{OUTDIR}/{sample}/kb_velo/matrix.ec'
    output:
        BCS = '{OUTDIR}/{sample}/kb_velo/spliced/output.barcodes.txt',
        GENES = '{OUTDIR}/{sample}/kb_velo/spliced/output.genes.txt',
        MAT = '{OUTDIR}/{sample}/kb_velo/spliced/output.mtx'
    params:
        MATDIR = directory('{OUTDIR}/{sample}/kb_velo/spliced')
    threads:
        1
    run:
        shell(
            f"""
            mkdir -p {params.MATDIR}

            {EXEC['BUSTOOLS']} count \
                --output {params.MATDIR}/ \
                --genemap {T2G_VELO_DICT[wildcards.sample]} \
                --ecmap {input.ECMAP} \
                --txnames {input.TRANSCRIPTS} \
                --genecounts \
                --umi-gene \
                --em \
                {input.SPLICED}
            """
        )

# unspliced matrix
rule bus2mat_velocity_unspliced:
    input:
        UNSPLICED = '{OUTDIR}/{sample}/kb_velo/unspliced.bus',
        TRANSCRIPTS = '{OUTDIR}/{sample}/kb_velo/transcripts.txt',
        ECMAP = '{OUTDIR}/{sample}/kb_velo/matrix.ec'
    output:
        BCS = '{OUTDIR}/{sample}/kb_velo/unspliced/output.barcodes.txt',
        GENES = '{OUTDIR}/{sample}/kb_velo/unspliced/output.genes.txt',
        MAT = '{OUTDIR}/{sample}/kb_velo/unspliced/output.mtx'
    params:
        MATDIR = directory('{OUTDIR}/{sample}/kb_velo/unspliced')
    threads:
        1
    run:
        shell(
            f"""
            mkdir -p {params.MATDIR}

            {EXEC['BUSTOOLS']} count \
                --output {params.MATDIR}/ \
                --genemap {T2G_VELO_DICT[wildcards.sample]} \
                --ecmap {input.ECMAP} \
                --txnames {input.TRANSCRIPTS} \
                --genecounts \
                --umi-gene \
                --em \
                {input.UNSPLICED}
            """
        )

# gzip the count matrix, etc.
rule compress_kb_outs_velocity:
    input:
        BCS_SPLICED = '{OUTDIR}/{sample}/kb_velo/spliced/output.barcodes.txt',
        GENES_SPLICED = '{OUTDIR}/{sample}/kb_velo/spliced/output.genes.txt',
        MAT_SPLICED = '{OUTDIR}/{sample}/kb_velo/spliced/output.mtx',
        BCS_UNSPLICED = '{OUTDIR}/{sample}/kb_velo/unspliced/output.barcodes.txt',
        GENES_UNSPLICED = '{OUTDIR}/{sample}/kb_velo/unspliced/output.genes.txt',
        MAT_UNSPLICED = '{OUTDIR}/{sample}/kb_velo/unspliced/output.mtx'
    output:
        BCS_SPLICED = '{OUTDIR}/{sample}/kb_velo/spliced/output.barcodes.txt.gz',
        GENES_SPLICED = '{OUTDIR}/{sample}/kb_velo/spliced/output.genes.txt.gz',
        MAT_SPLICED = '{OUTDIR}/{sample}/kb_velo/spliced/output.mtx.gz',
        BCS_UNSPLICED = '{OUTDIR}/{sample}/kb_velo/unspliced/output.barcodes.txt.gz',
        GENES_UNSPLICED = '{OUTDIR}/{sample}/kb_velo/unspliced/output.genes.txt.gz',
        MAT_UNSPLICED = '{OUTDIR}/{sample}/kb_velo/unspliced/output.mtx.gz'
    threads:
        config['CORES']        
    run:
        shell(
            f"""
            pigz -p{threads} \
                {input.BCS_SPLICED} {input.GENES_SPLICED} {input.MAT_SPLICED} \
                {input.BCS_UNSPLICED} {input.GENES_UNSPLICED} {input.MAT_UNSPLICED}
            """
        )