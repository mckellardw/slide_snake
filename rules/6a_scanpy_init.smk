# initialize & cache the **raw** counts as an anndata file for easier loading later
## Removes barcodes for which there are no molecules detected [`--remove_zero_features`]
rule cache_preQC_h5ad_STAR:
    input:
        BCS="{OUTDIR}/{SAMPLE}/STARsolo/{RECIPE}/Solo.out/{SOLO}/raw/barcodes.tsv.gz",
        FEATS="{OUTDIR}/{SAMPLE}/STARsolo/{RECIPE}/Solo.out/{SOLO}/raw/features.tsv.gz",
        MAT="{OUTDIR}/{SAMPLE}/STARsolo/{RECIPE}/Solo.out/{SOLO}/raw/{ALGO}.mtx.gz",
        BC_map="{OUTDIR}/{SAMPLE}/bc/whitelist_underscore.txt",
    output:
        H5AD="{OUTDIR}/{SAMPLE}/STARsolo/{RECIPE}/Solo.out/{SOLO}/raw/{ALGO}.h5ad",
    params:
        var_names="gene_symbols",  # scanpy.read_10x_mtx()
    threads: 1
    log:
        log="{OUTDIR}/{SAMPLE}/STARsolo/{RECIPE}/Solo.out/{SOLO}/raw/{ALGO}_cache.log",
    conda:
        f"{workflow.basedir}/envs/slsn_scanpy.yml"
    shell:
        """
        python scripts/py/cache_mtx_to_h5ad.py \
            --mat_in {input.MAT} \
            --feat_in {input.FEATS} \
            --bc_in {input.BCS} \
            --bb_map {input.BC_map} \
            --ad_out {output.H5AD} \
            --feat_col 1 \
            --remove_zero_features \
        1> {log.log}
        """


# initialize & cache the **raw** counts as an anndata file for easier loading later
## Removes barcodes for which there are no molecules detected [`--remove_zero_features`]
rule cache_preQC_h5ad_kb:
    input:
        BCS="{OUTDIR}/{SAMPLE}/kb/{RECIPE}/raw/output.barcodes.txt.gz",
        FEATS="{OUTDIR}/{SAMPLE}/kb/{RECIPE}/raw/output.genes.txt.gz",
        MAT="{OUTDIR}/{SAMPLE}/kb/{RECIPE}/raw/output.mtx.gz",
        BC_map="{OUTDIR}/{SAMPLE}/bc/whitelist_underscore.txt",
    output:
        H5AD="{OUTDIR}/{SAMPLE}/kb/{RECIPE}/raw/output.h5ad",
    log:
        log="{OUTDIR}/{SAMPLE}/kb/{RECIPE}/raw/cache.log",
    # params:
    #     var_names = "gene_symbols" # scanpy.read_10x_mtx()
    threads: 1
    conda:
        f"{workflow.basedir}/envs/slsn_scanpy.yml"
    shell:
        """
        python scripts/py/cache_mtx_to_h5ad.py \
            --mat_in {input.MAT} \
            --feat_in {input.FEATS} \
            --bc_in {input.BCS} \
            --bb_map {input.BC_map} \
            --ad_out {output.H5AD} \
            --feat_col 0 \
            --remove_zero_features \
        1> {log.log}
        """


rule cache_preQC_h5ad_kbpython:
    input:
        BCS="{OUTDIR}/{SAMPLE}/kbpython/{RECIPE}/counts_unfiltered/cells_x_genes.barcodes.txt.gz",
        FEATS="{OUTDIR}/{SAMPLE}/kbpython/{RECIPE}/counts_unfiltered/cells_x_genes.genes.txt.gz",
        MAT="{OUTDIR}/{SAMPLE}/kbpython/{RECIPE}/counts_unfiltered/cells_x_genes.mtx.gz",
        BC_map="{OUTDIR}/{SAMPLE}/bc/whitelist_underscore.txt",
    output:
        H5AD="{OUTDIR}/{SAMPLE}/kbpython/{RECIPE}/counts_unfiltered/output.h5ad",
    log:
        log="{OUTDIR}/{SAMPLE}/kbpython/{RECIPE}/counts_unfiltered/cache.log",
    # params:
    #     var_names = "gene_symbols" # scanpy.read_10x_mtx()
    threads: 1
    conda:
        f"{workflow.basedir}/envs/slsn_scanpy.yml"
    shell:
        """
        python scripts/py/cache_mtx_to_h5ad.py \
            --mat_in {input.MAT} \
            --feat_in {input.FEATS} \
            --bc_in {input.BCS} \
            --bb_map {input.BC_map} \
            --ad_out {output.H5AD} \
            --feat_col 0 \
            --remove_zero_features \
        1> {log.log}
        """


# initialize & cache the **raw** counts as an anndata file for easier loading later
## Removes barcodes for which there are no molecules detected [`--remove_zero_features`]
rule cache_preQC_h5ad_miRNA:
    input:
        BCS="{OUTDIR}/{SAMPLE}/miRNA/{RECIPE}/raw/barcodes.tsv.gz",
        FEATS="{OUTDIR}/{SAMPLE}/miRNA/{RECIPE}/{RECIPE}/raw/features.tsv.gz",
        MAT="{OUTDIR}/{SAMPLE}/miRNA/{RECIPE}/raw/matrix.mtx.gz",
        BC_map="{OUTDIR}/{SAMPLE}/bc/whitelist_underscore.txt",
    output:
        H5AD="{OUTDIR}/{SAMPLE}/miRNA/{RECIPE}/raw/output.h5ad",
    threads: 1
    conda:
        f"{workflow.basedir}/envs/slsn_scanpy.yml"
    shell:
        """
        python scripts/py/cache_mtx_to_h5ad.py \
            --mat_in {input.MAT} \
            --feat_in {input.FEATS} \
            --bc_in {input.BCS} \
            --bb_map {input.BC_map} \
            --ad_out {output.H5AD} \
            --feat_col 0 \
            --remove_zero_features
        """


# initialize & cache the **raw** counts as an anndata file for easier loading later
## Removes barcodes for which there are no molecules detected [`--remove_zero_features`]
rule cache_preQC_h5ad_piRNA:
    input:
        BCS="{OUTDIR}/{SAMPLE}/piRNA/{RECIPE}/raw/barcodes.tsv.gz",
        FEATS="{OUTDIR}/{SAMPLE}/piRNA/{RECIPE}/raw/features.tsv.gz",
        MAT="{OUTDIR}/{SAMPLE}/piRNA/{RECIPE}/raw/matrix.mtx.gz",
        BC_map="{OUTDIR}/{SAMPLE}/bc/whitelist_underscore.txt",
    output:
        H5AD="{OUTDIR}/{SAMPLE}/piRNA/{RECIPE}/raw/output.h5ad",
    threads: 1
    conda:
        f"{workflow.basedir}/envs/slsn_scanpy.yml"
    shell:
        """
        python scripts/py/cache_mtx_to_h5ad.py \
            --mat_in {input.MAT} \
            --feat_in {input.FEATS} \
            --bc_in {input.BCS} \
            --bb_map {input.BC_map} \
            --ad_out {output.H5AD} \
            --feat_col 0 \
            --remove_zero_features
        """


# initialize & cache the **raw** counts as an anndata file for easier loading later
## Removes barcodes for which there are no molecules detected [`--remove_zero_features`]
# rule ont_cache_preQC_h5ad_minimap2:
#     input:
#         MAT = "{OUTDIR}/{SAMPLE}/ont/minimap2/{RECIPE}/raw/umitools_counts.tsv.gz",
#         BC_map = lambda wildcards: BB_DICT[wildcards.SAMPLE] #TODO Adjust to match barcode handling schemas
#     output:
#         H5AD = "{OUTDIR}/{SAMPLE}/ont/minimap2/{RECIPE}/raw/output.h5ad"
#     threads:
#         1
#     conda:
#         f"{workflow.basedir}/envs/slsn_scanpy.yml"
#     shell:
#         """
#         python scripts/py/cache_umitools_to_h5ad.py \
#             --mat_in {input.MAT} \
#             --bb_map {input.BC_map}\
#             --ad_out {output.H5AD}\
#             --remove_zero_features
#         """


rule ont_cache_preQC_h5ad_minimap2:
    input:
        BCS="{OUTDIR}/{SAMPLE}/ont/minimap2/{RECIPE}/raw/barcodes.tsv.gz",
        FEATS="{OUTDIR}/{SAMPLE}/ont/minimap2/{RECIPE}/raw/features.tsv.gz",
        MAT="{OUTDIR}/{SAMPLE}/ont/minimap2/{RECIPE}/raw/matrix.mtx.gz",
        BC_map="{OUTDIR}/{SAMPLE}/bc/whitelist_underscore.txt",
    output:
        H5AD="{OUTDIR}/{SAMPLE}/ont/minimap2/{RECIPE}/raw/output.h5ad",
    log:
        log="{OUTDIR}/{SAMPLE}/ont/minimap2/{RECIPE}/raw/cache.log",
    threads: 1
    conda:
        f"{workflow.basedir}/envs/slsn_scanpy.yml"
    shell:
        """
        python scripts/py/cache_mtx_to_h5ad.py \
            --mat_in {input.MAT} \
            --feat_in {input.FEATS} \
            --bc_in {input.BCS} \
            --bb_map {input.BC_map} \
            --ad_out {output.H5AD} \
            --feat_col 0 \
            --remove_zero_features \
        1> {log.log}
        """
