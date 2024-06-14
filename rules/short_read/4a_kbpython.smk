#############################################
## kb-python pseudoalignment
##  source: https://github.com/pachterlab/kb_python
##  info: https://www.kallistobus.tools/  [currently very out of date]
#############################################
rule kbpython_std:
    input:
        FQS=lambda w: get_fqs(w, return_type="list", mode="ILMN"),
        BC="{OUTDIR}/{SAMPLE}/bc/whitelist.txt",
    output:
        BUS=temp("{OUTDIR}/{SAMPLE}/kbpython_std/{RECIPE}/output.unfiltered.bus"),
        # BUS_CORRECTED = temp('{OUTDIR}/{SAMPLE}/kbpython_std/{RECIPE}/output.corrected.bus'),
        # TRANSCRIPTS = '{OUTDIR}/{SAMPLE}/kbpython_std/{RECIPE}/transcripts.txt',
        ECMAP=temp("{OUTDIR}/{SAMPLE}/kbpython_std/{RECIPE}/matrix.ec"),
        # BCS="{OUTDIR}/{SAMPLE}/kbpython_std/{RECIPE}/counts_unfiltered/cells_x_genes.barcodes.txt",
        # FEATS="{OUTDIR}/{SAMPLE}/kbpython_std/{RECIPE}/counts_unfiltered/cells_x_genes.genes.txt",
        # MAT="{OUTDIR}/{SAMPLE}/kbpython_std/{RECIPE}/counts_unfiltered/cells_x_genes.mtx",        
        BCS="{OUTDIR}/{SAMPLE}/kbpython_std/{RECIPE}/counts_unfiltered/cellranger/barcodes.tsv",
        FEATS="{OUTDIR}/{SAMPLE}/kbpython_std/{RECIPE}/counts_unfiltered/cellranger/genes.tsv",
        MAT="{OUTDIR}/{SAMPLE}/kbpython_std/{RECIPE}/counts_unfiltered/cellranger/matrix.mtx",
    params:
        KB_IDX=lambda wildcards: IDX_DICT[wildcards.SAMPLE],
        KB_X=lambda wildcards: RECIPE_SHEET["kb.x"][wildcards.RECIPE],
        KB_EXTRA=lambda wildcards: RECIPE_SHEET["kb.extra"][wildcards.RECIPE],
        KB_T2G=lambda wildcards: T2G_DICT[wildcards.SAMPLE],
        N_READS_SUMMARY=1000000,  # number of reads to use for summary stats
        KB=EXEC["KB"],
        KALLISTO=EXEC["KALLISTO"],
        BUSTOOLS=EXEC["BUSTOOLS"],
    log:
        log="{OUTDIR}/{SAMPLE}/kbpython_std/{RECIPE}/kbpython_std.log",
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
            --workflow standard \
            -x '{params.KB_X}' \
            -w {input.BC} \
            --cellranger \
            --overwrite \
            -t {threads} {params.KB_EXTRA} \
            {input.FQS[0]} {input.FQS[1]} \
        2> {log.log}
        """
        # code to compute mean & std dev of read lengths...
        # --fragment-l $(zcat {input.FQS[1]} | head -n {params.N_READS_SUMMARY} | awk -f scripts/awk/fq_meanReadLength_int.awk) \
        # --fragment-s $(zcat {input.FQS[1]} | head -n {params.N_READS_SUMMARY} | awk -f scripts/awk/fq_stdDevReadLength_int.awk) \
        # meanLength=$(zcat {input.FQS[1]} | head -n {params.N_READS_SUMMARY} | awk -f scripts/awk/fq_meanReadLength.awk)
        # sdLength=$(zcat {input.FQS[1]} | head -n {params.N_READS_SUMMARY} | awk -f scripts/awk/fq_stdDevReadLength.awk)
        # echo "Using {params.N_READS_SUMMARY} for summary stats" > {log.log}
        # echo "Mean read length:               $meanLength" >> {log.log}
        # echo "Standard Deviation read length: $sdLength" >> {log.log}
        # -m {resources.MEM_GB} \

rule kbpython_std_remove_suffix:
    input:
        BCS="{OUTDIR}/{SAMPLE}/kbpython_std/{RECIPE}/counts_unfiltered/cellranger/barcodes.tsv",
    output:
        BCS="{OUTDIR}/{SAMPLE}/kbpython_std/{RECIPE}/counts_unfiltered/cellranger/barcodes_noSuffix.tsv"
    params:
        SUFFIX="-1"
    threads: config["CORES"]
    shell:
        """
        sed 's/{params.SUFFIX}//' {input.BCS} > {output.BCS}
        """

# # gzip the count matrix, etc.
rule kbpython_std_compress_outs:
    input:
        # BCS="{OUTDIR}/{SAMPLE}/kbpython_std/{RECIPE}/counts_unfiltered/cells_x_genes.barcodes.txt",
        # FEATS="{OUTDIR}/{SAMPLE}/kbpython_std/{RECIPE}/counts_unfiltered/cells_x_genes.genes.txt",
        # MAT="{OUTDIR}/{SAMPLE}/kbpython_std/{RECIPE}/counts_unfiltered/cells_x_genes.mtx",
        BCS="{OUTDIR}/{SAMPLE}/kbpython_std/{RECIPE}/counts_unfiltered/cellranger/barcodes.tsv",
        BCS2="{OUTDIR}/{SAMPLE}/kbpython_std/{RECIPE}/counts_unfiltered/cellranger/barcodes_noSuffix.tsv",
        FEATS="{OUTDIR}/{SAMPLE}/kbpython_std/{RECIPE}/counts_unfiltered/cellranger/genes.tsv",
        MAT="{OUTDIR}/{SAMPLE}/kbpython_std/{RECIPE}/counts_unfiltered/cellranger/matrix.mtx",
    output:
        # BCS="{OUTDIR}/{SAMPLE}/kbpython_std/{RECIPE}/counts_unfiltered/cells_x_genes.barcodes.txt.gz",
        # FEATS="{OUTDIR}/{SAMPLE}/kbpython_std/{RECIPE}/counts_unfiltered/cells_x_genes.genes.txt.gz",
        # MAT="{OUTDIR}/{SAMPLE}/kbpython_std/{RECIPE}/counts_unfiltered/cells_x_genes.mtx.gz",
        BCS="{OUTDIR}/{SAMPLE}/kbpython_std/{RECIPE}/counts_unfiltered/cellranger/barcodes.tsv.gz",
        BCS2="{OUTDIR}/{SAMPLE}/kbpython_std/{RECIPE}/counts_unfiltered/cellranger/barcodes_noSuffix.tsv.gz",
        FEATS="{OUTDIR}/{SAMPLE}/kbpython_std/{RECIPE}/counts_unfiltered/cellranger/genes.tsv.gz",
        MAT="{OUTDIR}/{SAMPLE}/kbpython_std/{RECIPE}/counts_unfiltered/cellranger/matrix.mtx.gz",
    threads: config["CORES"]
    shell:
        """
        pigz -p{threads} {input.BCS} {input.FEATS} {input.MAT} {input.BCS2}
        """
