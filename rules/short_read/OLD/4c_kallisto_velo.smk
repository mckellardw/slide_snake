#########################################################
## kallisto pseudoalignment for RNA velocity analysis
#########################################################
# Source: https://bustools.github.io/BUS_notebooks_R/velocity.html
# BUStools manual: https://bustools.github.io/manual


rule kallisto_align_velocity:
    input:
        R1_FQ="{OUTDIR}/{SAMPLE}/tmp/cut_R1.fq.gz",
        R2_FQ="{OUTDIR}/{SAMPLE}/tmp/cut_R2.fq.gz",
        R1_FQ_TWICE_CUT="{OUTDIR}/{SAMPLE}/tmp/twiceCut_R1.fq.gz",
        R2_FQ_TWICE_CUT="{OUTDIR}/{SAMPLE}/tmp/twiceCut_R2.fq.gz",
        R1_FQ_STAR_FILTERED="{OUTDIR}/{SAMPLE}/rRNA/STARsolo/final_filtered_R1.fq.gz",
        R2_FQ_STAR_FILTERED="{OUTDIR}/{SAMPLE}/rRNA/STARsolo/final_filtered_R2.fq.gz",
        R1_FQ_BWA_FILTERED="{OUTDIR}/{SAMPLE}/rRNA/bwa/final_filtered_R1.fq.gz",
        R2_FQ_BWA_FILTERED="{OUTDIR}/{SAMPLE}/rRNA/bwa/final_filtered_R2.fq.gz",
        BB="{OUTDIR}/{SAMPLE}/bc/whitelist.txt",
    output:
        BUS=temp("{OUTDIR}/{SAMPLE}/kb_velo/{RECIPE}/output.bus"),
        BUS_CORRECTED=temp("{OUTDIR}/{SAMPLE}/kb_velo/{RECIPE}/output.corrected.bus"),
        TRANSCRIPTS="{OUTDIR}/{SAMPLE}/kb_velo/{RECIPE}/transcripts.txt",
        ECMAP=temp("{OUTDIR}/{SAMPLE}/kb_velo/{RECIPE}/matrix.ec"),
    params:
        MEMLIMIT=config["MEMLIMIT_GB"],
    log:
        "{OUTDIR}/{SAMPLE}/kb_velo/{RECIPE}/kallisto_align.log",
    threads: config["CORES"]
    priority: 42
    run:
        # recipe = RECIPE_DICT[wildcards.SAMPLE]
        recipe = wildcards.RECIPE

        KB_IDX = IDX_DICT[wildcards.SAMPLE]
        KB_X = RECIPE_SHEET["kb.x"][recipe]

        # Select input reads based on alignment recipe
        if "rRNA.STAR" in recipe:  # Use trimmed & STAR-rRNA-filtered .fq's
            R1 = input.R1_FQ_STAR_FILTERED
            R2 = input.R2_FQ_STAR_FILTERED
        elif "rRNA.bwa" in recipe:  # TODO Use trimmed & bwa-rRNA-filtered .fq's
            R1 = input.R1_FQ_BWA_FILTERED
            R2 = input.R2_FQ_BWA_FILTERED
        elif "rRNA" not in recipe:  # just trimmed .fq's
            # R1 = input.R1_FQ
            # R2 = input.R2_FQ
            R1 = input.R1_FQ_TWICE_CUT
            R2 = input.R2_FQ_TWICE_CUT
        else:
            print("I just don't know what to do with myself...")

        shell(
            f"""
            bash scripts/bash/kb.sh \
                --outdir $(dirname {output.BUS}) \
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


# Split the .bus file for spliced/unspliced outputs
# TODO - fix hardcoded bits...
rule split_bus_velocity_spliced:
    input:
        BUS="{OUTDIR}/{SAMPLE}/kb_velo/{RECIPE}/output.corrected.bus",
        TRANSCRIPTS="{OUTDIR}/{SAMPLE}/kb_velo/{RECIPE}/transcripts.txt",
        ECMAP="{OUTDIR}/{SAMPLE}/kb_velo/{RECIPE}/matrix.ec",
    output:
        SPLICED="{OUTDIR}/{SAMPLE}/kb_velo/{RECIPE}/spliced.bus",
    log:
        log="{OUTDIR}/{SAMPLE}/kb_velo/{RECIPE}/split_bus_velocity_spliced.log",
    threads: 1
    run:
        KB_IDX = IDX_VELO_DICT[wildcards.SAMPLE]
        shell(
            f"""
            mkdir -p $(dirname {output.SPLICED})

            {EXEC['BUSTOOLS']} capture \
                --transcripts \
                --output {output.SPLICED} \
                --capture {KB_IDX}/cDNA.t2c \
                --ecmap {input.ECMAP} \
                --txnames {input.TRANSCRIPTS} \
                {input.BUS} \
            1> {log.log}
            """
        )


# Split the .bus file for spliced/unspliced outputs
# TODO - fix hardcoded bits...
rule split_bus_velocity_unspliced:
    input:
        BUS="{OUTDIR}/{SAMPLE}/kb_velo/{RECIPE}/output.corrected.bus",
        TRANSCRIPTS="{OUTDIR}/{SAMPLE}/kb_velo/{RECIPE}/transcripts.txt",
        ECMAP="{OUTDIR}/{SAMPLE}/kb_velo/{RECIPE}/matrix.ec",
    output:
        UNSPLICED="{OUTDIR}/{SAMPLE}/kb_velo/{RECIPE}/unspliced.bus",
    log:
        log="{OUTDIR}/{SAMPLE}/kb_velo/{RECIPE}/split_bus_velocity_unspliced.log",
    threads: 1
    run:
        KB_IDX = IDX_VELO_DICT[wildcards.SAMPLE]
        shell(
            f"""
            mkdir -p $(dirname {output.SPLICED})

            {EXEC['BUSTOOLS']} capture \
                --transcripts \
                --output {output.UNSPLICED} \
                --capture {KB_IDX}/introns.t2c \
                --ecmap {input.ECMAP} \
                --txnames {input.TRANSCRIPTS} \
                {input.BUS} \
            1> {log.log}
            """
        )


# build spliced matrix
rule bus2mat_velocity_spliced:
    input:
        SPLICED="{OUTDIR}/{SAMPLE}/kb_velo/{RECIPE}/spliced.bus",
        TRANSCRIPTS="{OUTDIR}/{SAMPLE}/kb_velo/{RECIPE}/transcripts.txt",
        ECMAP="{OUTDIR}/{SAMPLE}/kb_velo/{RECIPE}/matrix.ec",
    output:
        BCS="{OUTDIR}/{SAMPLE}/kb_velo/{RECIPE}/spliced/output.barcodes.txt",
        GENES="{OUTDIR}/{SAMPLE}/kb_velo/{RECIPE}/spliced/output.genes.txt",
        MAT="{OUTDIR}/{SAMPLE}/kb_velo/{RECIPE}/spliced/output.mtx",
    params:
        MATDIR=directory("{OUTDIR}/{SAMPLE}/kb_velo/{RECIPE}/spliced"),
    threads: 1
    run:
        shell(
            f"""
            mkdir -p $(dirname {output.MAT})

            {EXEC['BUSTOOLS']} count \
                --output $(dirname {output.MAT}) \
                --genemap {T2G_VELO_DICT[wildcards.SAMPLE]} \
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
        UNSPLICED="{OUTDIR}/{SAMPLE}/kb_velo/{RECIPE}/unspliced.bus",
        TRANSCRIPTS="{OUTDIR}/{SAMPLE}/kb_velo/{RECIPE}/transcripts.txt",
        ECMAP="{OUTDIR}/{SAMPLE}/kb_velo/{RECIPE}/matrix.ec",
    output:
        BCS="{OUTDIR}/{SAMPLE}/kb_velo/{RECIPE}/unspliced/output.barcodes.txt",
        GENES="{OUTDIR}/{SAMPLE}/kb_velo/{RECIPE}/unspliced/output.genes.txt",
        MAT="{OUTDIR}/{SAMPLE}/kb_velo/{RECIPE}/unspliced/output.mtx",
    params:
        MATDIR=directory("{OUTDIR}/{SAMPLE}/kb_velo/{RECIPE}/unspliced"),
    threads: 1
    run:
        shell(
            f"""
            mkdir -p $(dirname {output.MAT})

            {EXEC['BUSTOOLS']} count \
                --output $(dirname {output.MAT})/ \
                --genemap {T2G_VELO_DICT[wildcards.SAMPLE]} \
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
        BCS_SPLICED="{OUTDIR}/{SAMPLE}/kb_velo/{RECIPE}/spliced/output.barcodes.txt",
        GENES_SPLICED="{OUTDIR}/{SAMPLE}/kb_velo/{RECIPE}/spliced/output.genes.txt",
        MAT_SPLICED="{OUTDIR}/{SAMPLE}/kb_velo/{RECIPE}/spliced/output.mtx",
        BCS_UNSPLICED="{OUTDIR}/{SAMPLE}/kb_velo/{RECIPE}/unspliced/output.barcodes.txt",
        GENES_UNSPLICED="{OUTDIR}/{SAMPLE}/kb_velo/{RECIPE}/unspliced/output.genes.txt",
        MAT_UNSPLICED="{OUTDIR}/{SAMPLE}/kb_velo/{RECIPE}/unspliced/output.mtx",
    output:
        BCS_SPLICED="{OUTDIR}/{SAMPLE}/kb_velo/{RECIPE}/spliced/output.barcodes.txt.gz",
        GENES_SPLICED="{OUTDIR}/{SAMPLE}/kb_velo/{RECIPE}/spliced/output.genes.txt.gz",
        MAT_SPLICED="{OUTDIR}/{SAMPLE}/kb_velo/{RECIPE}/spliced/output.mtx.gz",
        BCS_UNSPLICED="{OUTDIR}/{SAMPLE}/kb_velo/{RECIPE}/unspliced/output.barcodes.txt.gz",
        GENES_UNSPLICED="{OUTDIR}/{SAMPLE}/kb_velo/{RECIPE}/unspliced/output.genes.txt.gz",
        MAT_UNSPLICED="{OUTDIR}/{SAMPLE}/kb_velo/{RECIPE}/unspliced/output.mtx.gz",
    threads: config["CORES"]
    run:
        shell(
            f"""
            {EXEC['PIGZ']} -p{threads} \
                {input.BCS_SPLICED} {input.GENES_SPLICED} {input.MAT_SPLICED} \
                {input.BCS_UNSPLICED} {input.GENES_UNSPLICED} {input.MAT_UNSPLICED}
            """
        )
