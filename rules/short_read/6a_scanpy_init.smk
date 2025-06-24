# initialize & cache the **raw** counts as an anndata file for easier loading later
## Removes barcodes for which there are no molecules detected [`--remove_zero_features`]


# STAR outputs
rule ilmn_6a_cache_h5ad_STAR:
    input:
        BCS="{OUTDIR}/{SAMPLE}/short_read/STARsolo/{RECIPE}/Solo.out/{SOLO}/raw/barcodes.tsv.gz",
        FEATS="{OUTDIR}/{SAMPLE}/short_read/STARsolo/{RECIPE}/Solo.out/{SOLO}/raw/features.tsv.gz",
        MAT="{OUTDIR}/{SAMPLE}/short_read/STARsolo/{RECIPE}/Solo.out/{SOLO}/raw/{ALGO}.mtx.gz",
        BC_map=lambda w: get_bc_map(w, mode="ILMN"),
    output:
        H5AD="{OUTDIR}/{SAMPLE}/short_read/STARsolo/{RECIPE}/Solo.out/{SOLO}/raw/{ALGO}.h5ad",
        QC_PLOTS="{OUTDIR}/{SAMPLE}/short_read/STARsolo/{RECIPE}/Solo.out/{SOLO}/raw/{ALGO}_qc_plots.png",
    params:
        var_names="gene_symbols",  # scanpy.read_10x_mtx()
    threads: 1
    log:
        log="{OUTDIR}/{SAMPLE}/short_read/STARsolo/{RECIPE}/Solo.out/{SOLO}/raw/{ALGO}_cache_h5ad.log",
        err="{OUTDIR}/{SAMPLE}/short_read/STARsolo/{RECIPE}/Solo.out/{SOLO}/raw/{ALGO}_cache_h5ad.err",
    conda:
        f"{workflow.basedir}/envs/scanpy.yml"
    shell:
        """
        python scripts/py/cache_mtx_to_h5ad.py \
            --mat_in {input.MAT} \
            --feat_in {input.FEATS} \
            --bc_in {input.BCS} \
            --bc_map {input.BC_map} \
            --ad_out {output.H5AD} \
            --feat_col 1 0 \
            --transpose \
            --remove_zero_features \
            --plot_qc \
            --qc_plot_file {output.QC_PLOTS} \
        1> {log.log} \
        2> {log.err}
        """


# kallisto/bustools outputs
rule ilmn_6a_cache_h5ad_kbpython_std:
    input:
        BCS="{OUTDIR}/{SAMPLE}/short_read/kbpython_std/{RECIPE}/counts_unfiltered/cellranger/barcodes_noSuffix.tsv.gz",
        FEATS="{OUTDIR}/{SAMPLE}/short_read/kbpython_std/{RECIPE}/counts_unfiltered/cellranger/genes.tsv.gz",
        MAT="{OUTDIR}/{SAMPLE}/short_read/kbpython_std/{RECIPE}/counts_unfiltered/cellranger/matrix.mtx.gz",
        BC_map="{OUTDIR}/{SAMPLE}/bc/map_underscore.txt",
    output:
        H5AD="{OUTDIR}/{SAMPLE}/short_read/kbpython_std/{RECIPE}/counts_unfiltered/output.h5ad",
        QC_PLOTS="{OUTDIR}/{SAMPLE}/short_read/kbpython_std/{RECIPE}/counts_unfiltered/qc_plots.png",
    log:
        log="{OUTDIR}/{SAMPLE}/short_read/kbpython_std/{RECIPE}/counts_unfiltered/cache_h5ad.log",
        err="{OUTDIR}/{SAMPLE}/short_read/kbpython_std/{RECIPE}/counts_unfiltered/cache_h5ad.err",
    threads: 1
    conda:
        f"{workflow.basedir}/envs/scanpy.yml"
    shell:
        """
        python scripts/py/cache_mtx_to_h5ad.py \
            --mat_in {input.MAT} \
            --feat_in {input.FEATS} \
            --bc_in {input.BCS} \
            --bc_map {input.BC_map} \
            --ad_out {output.H5AD} \
            --feat_col 1 0 \
            --transpose \
            --remove_zero_features \
            --plot_qc \
            --qc_plot_file {output.QC_PLOTS} \
        1> {log.log} \
        2> {log.err}
        """
