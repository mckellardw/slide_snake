#############################################
## kallisto pseudoalignment
#############################################
rule kallisto_align:
    input:
        # R1_FQ="{OUTDIR}/{SAMPLE}/tmp/cut_R1.fq.gz",
        # R2_FQ="{OUTDIR}/{SAMPLE}/tmp/cut_R2.fq.gz",
        # R1_FQ_TWICE_CUT="{OUTDIR}/{SAMPLE}/tmp/twiceCut_R1.fq.gz",
        # R2_FQ_TWICE_CUT="{OUTDIR}/{SAMPLE}/tmp/twiceCut_R2.fq.gz",
        # R1_FQ_STAR_FILTERED="{OUTDIR}/{SAMPLE}/rRNA/STARsolo/final_filtered_R1.fq.gz",
        # R2_FQ_STAR_FILTERED="{OUTDIR}/{SAMPLE}/rRNA/STARsolo/final_filtered_R2.fq.gz",
        # R1_FQ_BWA_FILTERED="{OUTDIR}/{SAMPLE}/rRNA/bwa/final_filtered_R1.fq.gz",
        # R2_FQ_BWA_FILTERED="{OUTDIR}/{SAMPLE}/rRNA/bwa/final_filtered_R2.fq.gz",
        FQS=lambda w: get_fqs(w, return_type="list", mode="ILMN"),
        BC="{OUTDIR}/{SAMPLE}/bc/whitelist.txt",
    output:
        BUS=temp("{OUTDIR}/{SAMPLE}/kb/{RECIPE}/output.bus"),
        BUS_CORRECTED=temp("{OUTDIR}/{SAMPLE}/kb/{RECIPE}/output.corrected.bus"),
        TRANSCRIPTS="{OUTDIR}/{SAMPLE}/kb/{RECIPE}/transcripts.txt",
        ECMAP=temp("{OUTDIR}/{SAMPLE}/kb/{RECIPE}/matrix.ec"),
    log:
        log="{OUTDIR}/{SAMPLE}/kb/{RECIPE}/kallisto_align.log",
    resources:
        threads=config["CORES"],
        mem_mb=config["MEMLIMIT_MB"],
        mem=config["MEMLIMIT_GB"],
    priority: 42
    run:
        recipe = wildcards.RECIPE

        KB_IDX = IDX_DICT[wildcards.SAMPLE]
        KB_X = RECIPE_SHEET["kb.x"][recipe]


        shell(
            f"""
            bash scripts/bash/kb.sh \
                --outdir $(dirname {output.BUS}) \
                --kb_idx {KB_IDX} \
                --whitelist {input.BC} \
                --chemistry {KB_X} \
                --log {log.log} \
                --threads {resources.threads} \
                --memlimit {params.mem} \
                --r1fq {input.FQS[0]} \
                --r2fq {input.FQS[1]}
            """
        )


rule bus2mat:
    input:
        BUS="{OUTDIR}/{SAMPLE}/kb/{RECIPE}/output.corrected.bus",
        TRANSCRIPTS="{OUTDIR}/{SAMPLE}/kb/{RECIPE}/transcripts.txt",
        ECMAP="{OUTDIR}/{SAMPLE}/kb/{RECIPE}/matrix.ec",
    output:
        BCS="{OUTDIR}/{SAMPLE}/kb/{RECIPE}/raw/output.barcodes.txt",
        GENES="{OUTDIR}/{SAMPLE}/kb/{RECIPE}/raw/output.genes.txt",
        MAT="{OUTDIR}/{SAMPLE}/kb/{RECIPE}/raw/output.mtx",
    log:
        log="{OUTDIR}/{SAMPLE}/kb/{RECIPE}/raw/bustools_count.log",
    threads: 1
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
            2> {log.log}
            """
        )


# gzip the count matrix, etc.
rule compress_kb_outs:
    input:
        BCS="{OUTDIR}/{SAMPLE}/kb/{RECIPE}/raw/output.barcodes.txt",
        GENES="{OUTDIR}/{SAMPLE}/kb/{RECIPE}/raw/output.genes.txt",
        MAT="{OUTDIR}/{SAMPLE}/kb/{RECIPE}/raw/output.mtx",
    output:
        BCS="{OUTDIR}/{SAMPLE}/kb/{RECIPE}/raw/output.barcodes.txt.gz",
        GENES="{OUTDIR}/{SAMPLE}/kb/{RECIPE}/raw/output.genes.txt.gz",
        MAT="{OUTDIR}/{SAMPLE}/kb/{RECIPE}/raw/output.mtx.gz",
    threads: config["CORES"]
    run:
        shell(
            f"""
            {EXEC['PIGZ']} -p{resources.threads} {input.BCS} {input.GENES} {input.MAT}
            """
        )
