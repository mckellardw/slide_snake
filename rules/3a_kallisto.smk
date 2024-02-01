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
        # recipe = RECIPE_DICT[wildcards.SAMPLE]
        recipe = wildcards.RECIPE
        
        KB_IDX = IDX_DICT[wildcards.SAMPLE]        
        KB_X = RECIPE_SHEET["kb.x"][recipe]
        
        # Select input reads based on alignment recipe
        if "rRNA.STAR" in recipe: # Use trimmed & STAR-rRNA-filtered .fq's
            R1 = input.R1_FQ_FILTERED
            R2 = input.R2_FQ_FILTERED
        elif "rRNA.bwa" in recipe: #TODO Use trimmed & bwa-rRNA-filtered .fq's
            print("TODO")
            # R1 = input.R1_FQ_FILTERED
            # R2 = input.R2_FQ_FILTERED
        elif "rRNA" not in recipe: # just trimmed .fq's
            R1 = input.R1_FQ
            R2 = input.R2_FQ
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

rule bus2mat:
    input:
        BUS         = '{OUTDIR}/{SAMPLE}/kb/{RECIPE}/output.corrected.bus',
        TRANSCRIPTS = '{OUTDIR}/{SAMPLE}/kb/{RECIPE}/transcripts.txt',
        ECMAP       = '{OUTDIR}/{SAMPLE}/kb/{RECIPE}/matrix.ec'
    output:
        BCS   = '{OUTDIR}/{SAMPLE}/kb/{RECIPE}/raw/output.barcodes.txt',
        GENES = '{OUTDIR}/{SAMPLE}/kb/{RECIPE}/raw/output.genes.txt',
        MAT   = '{OUTDIR}/{SAMPLE}/kb/{RECIPE}/raw/output.mtx'
    params:
    log:
        '{OUTDIR}/{SAMPLE}/kb/{RECIPE}/raw/bustools_count.log'
    threads:
        1
    run:
        KB_T2G = T2G_DICT[wildcards.SAMPLE]

        shell(
            f"""
            mkdir -p $(dirname {output.MAT})

            {EXEC['BUSTOOLS']} count \
                --output $(dirname {output.MAT})/ \
                --genemap {KB_T2G} \
                --ecmap {input.ECMAP} \
                --txnames {input.TRANSCRIPTS} \
                --genecounts \
                --umi-gene \
                --em \
                {input.BUS} \
            2> {log}
            """
        )

#TODO- use `--downsample` to generate count mats with different downsampling rates
# rule kb_rarefaction:
#     input:
#         BUS         = '{OUTDIR}/{SAMPLE}/kb/{RECIPE}/output.corrected.bus',
#         TRANSCRIPTS = '{OUTDIR}/{SAMPLE}/kb/{RECIPE}/transcripts.txt',
#         ECMAP       = '{OUTDIR}/{SAMPLE}/kb/{RECIPE}/matrix.ec'
#     output:
#         BCS   = '{OUTDIR}/{SAMPLE}/kb/{RECIPE}/raw/output.barcodes.txt',
#         GENES = '{OUTDIR}/{SAMPLE}/kb/{RECIPE}/raw/output.genes.txt',
#         MAT   = '{OUTDIR}/{SAMPLE}/kb/{RECIPE}/raw/output.mtx'
#     params:
#         DOWN_RATES = [0.1, 0.2, 0.5, 0.75]
#     log:
#         '{OUTDIR}/{SAMPLE}/kb/{RECIPE}/raw/rarefaction.log'
#     threads:
#         1
#     run:
#         KB_T2G = T2G_DICT[wildcards.SAMPLE]

#         shell(
#             f"""
#             mkdir -p {params.MATDIR}

#             {EXEC['BUSTOOLS']} count \
#                 --output $(dirname {output.MAT}) \
#                 --genemap {KB_T2G} \
#                 --ecmap {input.ECMAP} \
#                 --txnames {input.TRANSCRIPTS} \
#                 --genecounts \
#                 --umi-gene \
#                 --em \
#                 --downsample {DOWN_RATE}
#                 {input.BUS} \
#             2> {log}
#             """
#         )


# gzip the count matrix, etc.
rule compress_kb_outs:
    input:
        BCS   = '{OUTDIR}/{SAMPLE}/kb/{RECIPE}/raw/output.barcodes.txt',
        GENES = '{OUTDIR}/{SAMPLE}/kb/{RECIPE}/raw/output.genes.txt',
        MAT   = '{OUTDIR}/{SAMPLE}/kb/{RECIPE}/raw/output.mtx'
    output:
        BCS   = '{OUTDIR}/{SAMPLE}/kb/{RECIPE}/raw/output.barcodes.txt.gz',
        GENES = '{OUTDIR}/{SAMPLE}/kb/{RECIPE}/raw/output.genes.txt.gz',
        MAT   = '{OUTDIR}/{SAMPLE}/kb/{RECIPE}/raw/output.mtx.gz'
    params:
    threads:
        config['CORES']        
    run:
        shell(
            f"""
            {EXEC['PIGZ']} -p{threads} {input.BCS} {input.GENES} {input.MAT}
            """
        )