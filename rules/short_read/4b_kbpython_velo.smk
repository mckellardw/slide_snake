#############################################
## kb-python pseudoalignment
##  source: https://github.com/pachterlab/kb_python
##  info: https://www.kallistobus.tools/  [currently very out of date]
#############################################
rule kbpython_nac:
    input:
        FQS=lambda w: get_fqs(w, return_type="list", mode="ILMN"),
        BC="{OUTDIR}/{SAMPLE}/bc/whitelist.txt",
    output:
        BUS=temp("{OUTDIR}/{SAMPLE}/kbpython_nac/{RECIPE}/output.unfiltered.bus"),
        # BUS_CORRECTED = temp('{OUTDIR}/{SAMPLE}/kbpython_nac/{RECIPE}/output.corrected.bus'),
        # TRANSCRIPTS = '{OUTDIR}/{SAMPLE}/kbpython_nac/{RECIPE}/transcripts.txt',
        ECMAP=temp("{OUTDIR}/{SAMPLE}/kbpython_nac/{RECIPE}/matrix.ec"),
        BCS="{OUTDIR}/{SAMPLE}/kbpython_nac/{RECIPE}/counts_unfiltered/cells_x_genes.barcodes.txt",
        GENES="{OUTDIR}/{SAMPLE}/kbpython_nac/{RECIPE}/counts_unfiltered/cells_x_genes.genes.txt",
        MAT="{OUTDIR}/{SAMPLE}/kbpython_nac/{RECIPE}/counts_unfiltered/cells_x_genes.mtx",
    params:
        KB_X=lambda wildcards: RECIPE_SHEET["kb.x"][wildcards.RECIPE],
        KB_EXTRA=lambda wildcards: RECIPE_SHEET["kb.extra"][wildcards.RECIPE],
        KB_IDX=lambda wildcards: IDX_VELO_DICT[wildcards.SAMPLE],
        KB_T2G=lambda wildcards: T2G_DICT[wildcards.SAMPLE],
        KALLISTO=EXEC["KALLISTO"],
        KB=EXEC["KB"],
    log:
        log="{OUTDIR}/{SAMPLE}/kbpython_nac/{RECIPE}/kbpython_nac.log",
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
            -o $(dirname {output.BUS}) \
            --strand forward \
            --workflow nac \
            -x {params.KB_X} \
            -w {input.BC} \
            --overwrite \
            -t {threads} {params.KB_EXTRA} \
            {params.FQS[0]} {params.FQS[1]} \
        2> {log.log}
        """
        # -m {resources.MEM_GB} \
        # --report \


# # gzip the count matrix, etc.
rule kbpython_nac_compress_outs:
    input:
        BCS="{OUTDIR}/{SAMPLE}/kbpython_nac/{RECIPE}/counts_unfiltered/cells_x_genes.barcodes.txt",
        GENES="{OUTDIR}/{SAMPLE}/kbpython_nac/{RECIPE}/counts_unfiltered/cells_x_genes.genes.txt",
        MAT="{OUTDIR}/{SAMPLE}/kbpython_nac/{RECIPE}/counts_unfiltered/cells_x_genes.mtx",
    output:
        BCS="{OUTDIR}/{SAMPLE}/kbpython_nac/{RECIPE}/counts_unfiltered/cells_x_genes.barcodes.txt.gz",
        GENES="{OUTDIR}/{SAMPLE}/kbpython_nac/{RECIPE}/counts_unfiltered/cells_x_genes.genes.txt.gz",
        MAT="{OUTDIR}/{SAMPLE}/kbpython_nac/{RECIPE}/counts_unfiltered/cells_x_genes.mtx.gz",
    threads: config["CORES"]
    run:
        shell(
            f"""
            {EXEC['PIGZ']} -p{threads} {input.BCS} {input.GENES} {input.MAT}
            """
        )
