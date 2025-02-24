# initialize & cache the **raw** counts as a Seurat object for easier loading later
## Removes barcodes for which there are no molecules detected


# STAR outputs
rule ilmn_6a_cache_seurat_STAR:
    input:
        BCS="{OUTDIR}/{SAMPLE}/short_read/STARsolo/{RECIPE}/Solo.out/{SOLO}/raw/barcodes.tsv.gz",
        FEATS="{OUTDIR}/{SAMPLE}/short_read/STARsolo/{RECIPE}/Solo.out/{SOLO}/raw/features.tsv.gz",
        MAT="{OUTDIR}/{SAMPLE}/short_read/STARsolo/{RECIPE}/Solo.out/{SOLO}/raw/{ALGO}.mtx.gz",
        BC_map=lambda w: get_bc_map(w, mode="ILMN"),
    output:
        SEURAT="{OUTDIR}/{SAMPLE}/short_read/STARsolo/{RECIPE}/Solo.out/{SOLO}/raw/{ALGO}.rds",
    params:
        var_names="gene_symbols",  # Seurat::CreateSeuratObject()
    threads: 1
    log:
        log="{OUTDIR}/{SAMPLE}/short_read/STARsolo/{RECIPE}/Solo.out/{SOLO}/raw/{ALGO}_cache_rds.log",
        err="{OUTDIR}/{SAMPLE}/short_read/STARsolo/{RECIPE}/Solo.out/{SOLO}/raw/{ALGO}_cache_rds.err",
    conda:
        f"{workflow.basedir}/envs/seurat.yml"
    shell:
        """
        Rscript scripts/R/cache_mtx_to_seurat.R \
            --mat_in {input.MAT} \
            --feat_in {input.FEATS} \
            --bc_in {input.BCS} \
            --bc_map {input.BC_map} \
            --seurat_out {output.SEURAT} \
        1> {log.log} \
        2> {log.err}
        """


# kallisto/bustools outputs
rule ilmn_6a_cache_seurat_kbpython_std:
    input:
        BCS="{OUTDIR}/{SAMPLE}/short_read/kbpython_std/{RECIPE}/counts_unfiltered/cellranger/barcodes_noSuffix.tsv.gz",
        FEATS="{OUTDIR}/{SAMPLE}/short_read/kbpython_std/{RECIPE}/counts_unfiltered/cellranger/genes.tsv.gz",
        MAT="{OUTDIR}/{SAMPLE}/short_read/kbpython_std/{RECIPE}/counts_unfiltered/cellranger/matrix.mtx.gz",
        BC_map="{OUTDIR}/{SAMPLE}/bc/map_underscore.txt",
    output:
        SEURAT="{OUTDIR}/{SAMPLE}/short_read/kbpython_std/{RECIPE}/counts_unfiltered/output.rds",
    log:
        log="{OUTDIR}/{SAMPLE}/short_read/kbpython_std/{RECIPE}/counts_unfiltered/cache_rds.log",
        err="{OUTDIR}/{SAMPLE}/short_read/kbpython_std/{RECIPE}/counts_unfiltered/cache_rds.err",
    threads: 1
    conda:
        f"{workflow.basedir}/envs/seurat.yml"
    shell:
        """
        Rscript scripts/R/cache_mtx_to_seurat.R \
            --mat_in {input.MAT} \
            --feat_in {input.FEATS} \
            --bc_in {input.BCS} \
            --bc_map {input.BC_map} \
            --seurat_out {output.SEURAT} \
        1> {log.log} \
        2> {log.err}
        """
