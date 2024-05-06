#############################################
## kb-python pseudoalignment
##  source: https://github.com/pachterlab/kb_python
##  info: https://www.kallistobus.tools/  [currently very out of date]
#############################################
rule kbpython_velo:
    input:
        R1_FQ="{OUTDIR}/{SAMPLE}/tmp/cut_R1.fq.gz",
        R2_FQ="{OUTDIR}/{SAMPLE}/tmp/cut_R2.fq.gz",
        R1_FQ_TWICE_CUT="{OUTDIR}/{SAMPLE}/tmp/twiceCut_R1.fq.gz",
        R2_FQ_TWICE_CUT="{OUTDIR}/{SAMPLE}/tmp/twiceCut_R2.fq.gz",
        R1_FQ_STAR_FILTERED="{OUTDIR}/{SAMPLE}/rRNA/STARsolo/final_filtered_R1.fq.gz",
        R2_FQ_STAR_FILTERED="{OUTDIR}/{SAMPLE}/rRNA/STARsolo/final_filtered_R2.fq.gz",
        R1_FQ_BWA_FILTERED="{OUTDIR}/{SAMPLE}/rRNA/bwa/final_filtered_R1.fq.gz",
        R2_FQ_BWA_FILTERED="{OUTDIR}/{SAMPLE}/rRNA/bwa/final_filtered_R2.fq.gz",
        BC="{OUTDIR}/{SAMPLE}/bb/whitelist.txt",
    output:
        BUS=temp("{OUTDIR}/{SAMPLE}/kbpython/{RECIPE}/output.unfiltered.bus"),
        # BUS_CORRECTED = temp('{OUTDIR}/{SAMPLE}/kbpython/{RECIPE}/output.corrected.bus'),
        # TRANSCRIPTS = '{OUTDIR}/{SAMPLE}/kbpython/{RECIPE}/transcripts.txt',
        ECMAP=temp("{OUTDIR}/{SAMPLE}/kbpython/{RECIPE}/matrix.ec"),
        BCS="{OUTDIR}/{SAMPLE}/kbpython/{RECIPE}/counts_unfiltered/cells_x_genes.barcodes.txt",
        GENES="{OUTDIR}/{SAMPLE}/kbpython/{RECIPE}/counts_unfiltered/cells_x_genes.genes.txt",
        MAT="{OUTDIR}/{SAMPLE}/kbpython/{RECIPE}/counts_unfiltered/cells_x_genes.mtx",
    params:
        KB_X=lambda wildcards: RECIPE_SHEET["kb.x"][wildcards.RECIPE],
        KB_IDX=lambda wildcards: IDX_VELO_DICT[wildcards.SAMPLE],
        KB_T2G=lambda wildcards: T2G_DICT[wildcards.SAMPLE],
        FQS=lambda w: get_fqs(w),
        KALLISTO=EXEC["KALLISTO"],
        KB=EXEC["KB"],
    log:
        log="{OUTDIR}/{SAMPLE}/kbpython/{RECIPE}/kbpython_standard.log",
    threads: config["CORES"]
    resources:
        MEM_GB=config["MEMLIMIT_GB"],
    priority: 42
    conda:
        f"{workflow.basedir}/envs/kb.yml"
    shell:
        """
        mkdir -p $(dirname {output.BUS})

        {params.KB} count \
            -i {params.KB_IDX} \
            -g {params.KB_T2G} \
            --kallisto {params.KALLISTO} \
            --bustools {params.BUSTOOLS} \
            -o $(dirname {output.BUS}) \
            --strand forward \
            -mm \
            --workflow nac \
            -x {params.KB_X} \
            -w {input.BC} \
            -t {threads} \
            {params.FQS[0]} {params.FQS[1]} \
        2> {log.log}
        """
        # -m {resources.MEM_GB} \
        # --report \


# # gzip the count matrix, etc.
rule compress_kbpython_outs:
    input:
        BCS="{OUTDIR}/{SAMPLE}/kbpython/{RECIPE}/counts_unfiltered/cells_x_genes.barcodes.txt",
        GENES="{OUTDIR}/{SAMPLE}/kbpython/{RECIPE}/counts_unfiltered/cells_x_genes.genes.txt",
        MAT="{OUTDIR}/{SAMPLE}/kbpython/{RECIPE}/counts_unfiltered/cells_x_genes.mtx",
    output:
        BCS="{OUTDIR}/{SAMPLE}/kbpython/{RECIPE}/counts_unfiltered/cells_x_genes.barcodes.txt.gz",
        GENES="{OUTDIR}/{SAMPLE}/kbpython/{RECIPE}/counts_unfiltered/cells_x_genes.genes.txt.gz",
        MAT="{OUTDIR}/{SAMPLE}/kbpython/{RECIPE}/counts_unfiltered/cells_x_genes.mtx.gz",
    threads: config["CORES"]
    run:
        shell(
            f"""
            {EXEC['PIGZ']} -p{threads} {input.BCS} {input.GENES} {input.MAT}
            """
        )
