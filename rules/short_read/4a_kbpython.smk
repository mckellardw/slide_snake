#############################################
## kb-python pseudoalignment
##  source: https://github.com/pachterlab/kb_python
##  info: https://www.kallistobus.tools/  [currently very out of date]
#############################################
rule kbpython_standard:
    input:
        FQS=lambda w: get_fqs(w, return_type="list", mode="ILMN"),
        BC="{OUTDIR}/{SAMPLE}/bc/whitelist.txt",
    output:
        BUS=temp("{OUTDIR}/{SAMPLE}/kbpython/{RECIPE}/output.unfiltered.bus"),
        # BUS_CORRECTED = temp('{OUTDIR}/{SAMPLE}/kbpython/{RECIPE}/output.corrected.bus'),
        # TRANSCRIPTS = '{OUTDIR}/{SAMPLE}/kbpython/{RECIPE}/transcripts.txt',
        ECMAP=temp("{OUTDIR}/{SAMPLE}/kbpython/{RECIPE}/matrix.ec"),
        BCS="{OUTDIR}/{SAMPLE}/kbpython/{RECIPE}/counts_unfiltered/cells_x_genes.barcodes.txt",
        GENES="{OUTDIR}/{SAMPLE}/kbpython/{RECIPE}/counts_unfiltered/cells_x_genes.genes.txt",
        MAT="{OUTDIR}/{SAMPLE}/kbpython/{RECIPE}/counts_unfiltered/cells_x_genes.mtx",
    params:
        KB_IDX=lambda wildcards: IDX_DICT[wildcards.SAMPLE],
        KB_X=lambda wildcards: RECIPE_SHEET["kb.x"][wildcards.RECIPE],
        KB_T2G=lambda wildcards: T2G_DICT[wildcards.SAMPLE],
        KB=EXEC["KB"],
        KALLISTO=EXEC["KALLISTO"],
        BUSTOOLS=EXEC["BUSTOOLS"],
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

        kb count \
            -i {params.KB_IDX} \
            -g {params.KB_T2G} \
            -o $(dirname {output.BUS}) \
            --strand forward \
            -mm \
            --workflow standard \
            -x {params.KB_X} \
            -w {input.BC} \
            -t {threads} \
            {input.FQS[0]} {input.FQS[1]} \
        2> {log.log}
        """
        # {params.KB}
            # --kallisto {params.KALLISTO} \
            # --bustools {params.BUSTOOLS} \


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
