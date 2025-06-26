#############################################
## kb-python pseudoalignment
##  source: https://github.com/pachterlab/kb_python
##  info: https://www.kallistobus.tools/  [currently very out of date]
#############################################
rule ilmn_4a_kbpython_std:
    input:
        FQS=lambda w: get_fqs(w, return_type="list", mode="ILMN"),
        BC="{OUTDIR}/{SAMPLE}/bc/whitelist.txt",
    output:
        BUS=temp("{OUTDIR}/{SAMPLE}/short_read/kbpython_std/{RECIPE}/output.bus"),
        ECMAP=temp("{OUTDIR}/{SAMPLE}/short_read/kbpython_std/{RECIPE}/matrix.ec"),
        BCS="{OUTDIR}/{SAMPLE}/short_read/kbpython_std/{RECIPE}/counts_unfiltered/cellranger/barcodes.tsv",
        FEATS="{OUTDIR}/{SAMPLE}/short_read/kbpython_std/{RECIPE}/counts_unfiltered/cellranger/genes.tsv",
        MAT="{OUTDIR}/{SAMPLE}/short_read/kbpython_std/{RECIPE}/counts_unfiltered/cellranger/matrix.mtx",
    params:
        KB_IDX=lambda w: get_kallisto_ref(w, mode="idx"),
        KB_T2G=lambda w: get_kallisto_ref(w, mode="t2g"),
        KB_X=lambda w: RECIPE_SHEET["kb_x"][w.RECIPE],
        KB_EXTRA=lambda w: RECIPE_SHEET["kb_extra"][w.RECIPE],
        N_READS_SUMMARY=1000000,  # number of reads to use for summary stats
    log:
        log="{OUTDIR}/{SAMPLE}/short_read/kbpython_std/{RECIPE}/kbpython_std.log",
        err="{OUTDIR}/{SAMPLE}/short_read/kbpython_std/{RECIPE}/kbpython_std.err",
    resources:
        mem=config["MEMLIMIT_GB"],
    threads: config["CORES"]
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
            -m 32G \
            -t {threads} {params.KB_EXTRA} \
            {input.FQS[0]} {input.FQS[1]} \
        1> {log.log} \
        2> {log.err}
        """


rule ilmn_4a_kbpython_std_remove_suffix:
    input:
        BCS="{OUTDIR}/{SAMPLE}/short_read/kbpython_std/{RECIPE}/counts_unfiltered/cellranger/barcodes.tsv",
    output:
        BCS="{OUTDIR}/{SAMPLE}/short_read/kbpython_std/{RECIPE}/counts_unfiltered/cellranger/barcodes_noSuffix.tsv",
    params:
        SUFFIX="-1",
    threads: 1
    shell:
        """
        sed 's/{params.SUFFIX}//' {input.BCS} > {output.BCS}
        """


# # gzip the count matrix, etc.
rule ilmn_4a_kbpython_std_compress_outs:
    input:
        BCS="{OUTDIR}/{SAMPLE}/short_read/kbpython_std/{RECIPE}/counts_unfiltered/cellranger/barcodes.tsv",
        BCS2="{OUTDIR}/{SAMPLE}/short_read/kbpython_std/{RECIPE}/counts_unfiltered/cellranger/barcodes_noSuffix.tsv",
        FEATS="{OUTDIR}/{SAMPLE}/short_read/kbpython_std/{RECIPE}/counts_unfiltered/cellranger/genes.tsv",
        MAT="{OUTDIR}/{SAMPLE}/short_read/kbpython_std/{RECIPE}/counts_unfiltered/cellranger/matrix.mtx",
    output:
        BCS="{OUTDIR}/{SAMPLE}/short_read/kbpython_std/{RECIPE}/counts_unfiltered/cellranger/barcodes.tsv.gz",
        BCS2="{OUTDIR}/{SAMPLE}/short_read/kbpython_std/{RECIPE}/counts_unfiltered/cellranger/barcodes_noSuffix.tsv.gz",
        FEATS="{OUTDIR}/{SAMPLE}/short_read/kbpython_std/{RECIPE}/counts_unfiltered/cellranger/genes.tsv.gz",
        MAT="{OUTDIR}/{SAMPLE}/short_read/kbpython_std/{RECIPE}/counts_unfiltered/cellranger/matrix.mtx.gz",
    threads: config["CORES"]
    shell:
        """
        pigz -p{threads} {input.BCS} {input.FEATS} {input.MAT} {input.BCS2}
        """
