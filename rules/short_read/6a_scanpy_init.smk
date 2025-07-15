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
        GTF=lambda w: SAMPLE_SHEET["genes_gtf"][w.SAMPLE],
        GTF_FEATURE_TYPE="gene",  # feature type in gtf to use 
        GTF_ID="gene_id",  # gtf attribute used to match var_names in adata
        FEAT_COL=0,  # column in features.tsv to use as var_names
    threads: 1
    log:
        log="{OUTDIR}/{SAMPLE}/short_read/STARsolo/{RECIPE}/Solo.out/{SOLO}/raw/logs/{ALGO}_cache_h5ad.log",
        err="{OUTDIR}/{SAMPLE}/short_read/STARsolo/{RECIPE}/Solo.out/{SOLO}/raw/logs/{ALGO}_cache_h5ad.err",
    conda:
        f"{workflow.basedir}/envs/scanpy.yml"
    shell:
        """
        mkdir -p $(dirname {log.log})
        python scripts/py/cache_mtx_to_h5ad.py \
            --mat_in {input.MAT} \
            --feat_in {input.FEATS} \
            --bc_in {input.BCS} \
            --bc_map {input.BC_map} \
            --ad_out {output.H5AD} \
            --feat_col {params.FEAT_COL} \
            --remove_zero_features \
            --plot_qc \
            --qc_plot_file {output.QC_PLOTS} \
            --gtf_file {params.GTF} \
            --gtf_feature_type {params.GTF_FEATURE_TYPE} \
            --gtf_id {params.GTF_ID} \
            --transpose \
        1> {log.log} \
        2> {log.err}
        """


# kallisto/bustools outputs
rule ilmn_6a_cache_h5ad_kbpython_std:
    input:
        BCS="{OUTDIR}/{SAMPLE}/short_read/kbpython_std/{RECIPE}/counts_unfiltered/cellranger/barcodes_noSuffix.tsv.gz",
        FEATS="{OUTDIR}/{SAMPLE}/short_read/kbpython_std/{RECIPE}/counts_unfiltered/cellranger/genes.tsv.gz",
        MAT="{OUTDIR}/{SAMPLE}/short_read/kbpython_std/{RECIPE}/counts_unfiltered/cellranger/matrix.mtx.gz",
        BC_map="{OUTDIR}/{SAMPLE}/bc/map.txt",
    output:
        H5AD="{OUTDIR}/{SAMPLE}/short_read/kbpython_std/{RECIPE}/counts_unfiltered/output.h5ad",
        QC_PLOTS="{OUTDIR}/{SAMPLE}/short_read/kbpython_std/{RECIPE}/counts_unfiltered/qc_plots.png",
    params:
        GTF=lambda w: SAMPLE_SHEET["genes_gtf"][w.SAMPLE],
        GTF_FEATURE_TYPE="gene",  # feature type in gtf to use 
        GTF_ID="gene_id",  # gtf attribute used to match var_names in adata
        FEAT_COL=0,  # column in features.tsv to use as var_names
    log:
        log="{OUTDIR}/{SAMPLE}/short_read/kbpython_std/{RECIPE}/counts_unfiltered/logs/cache_h5ad.log",
        err="{OUTDIR}/{SAMPLE}/short_read/kbpython_std/{RECIPE}/counts_unfiltered/logs/cache_h5ad.err",
    threads: 1
    conda:
        f"{workflow.basedir}/envs/scanpy.yml"
    shell:
        """
        mkdir -p $(dirname {log.log})
        python scripts/py/cache_mtx_to_h5ad.py \
            --mat_in {input.MAT} \
            --feat_in {input.FEATS} \
            --bc_in {input.BCS} \
            --bc_map {input.BC_map} \
            --ad_out {output.H5AD} \
            --feat_col {params.FEAT_COL} \
            --remove_zero_features \
            --plot_qc \
            --qc_plot_file {output.QC_PLOTS} \
            --gtf_file {params.GTF} \
            --gtf_feature_type {params.GTF_FEATURE_TYPE} \
            --gtf_id {params.GTF_ID} \
            --transpose \
        1> {log.log} \
        2> {log.err}
        """
