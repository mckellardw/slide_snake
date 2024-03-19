# initialize & cache the **raw** counts as an anndata file for easier loading later
## Removes barcodes for which there are no molecules detected [`--remove_zero_features`]
rule cache_preQC_h5ad_STAR:
    input:
        BCS = '{OUTDIR}/{SAMPLE}/STARsolo/{RECIPE}/Solo.out/{SOLO}/raw/barcodes.tsv.gz',
        GENES = '{OUTDIR}/{SAMPLE}/STARsolo/{RECIPE}/Solo.out/{SOLO}/raw/features.tsv.gz',
        MAT = '{OUTDIR}/{SAMPLE}/STARsolo/{RECIPE}/Solo.out/{SOLO}/raw/{ALGO}.mtx.gz',
        BB_map = lambda wildcards: BB_DICT[wildcards.SAMPLE]
    output:
        H5AD = '{OUTDIR}/{SAMPLE}/STARsolo/{RECIPE}/Solo.out/{SOLO}/raw/{ALGO}.h5ad'
    params:
        var_names = "gene_symbols" # scanpy.read_10x_mtx()
    threads:
        1
    run:
        shell(
            f"""
            python scripts/py/cache_mtx_to_h5ad.py \
                --mat_in {input.MAT} \
                --feat_in {input.GENES} \
                --bc_in {input.BCS} \
                --bb_map {input.BB_map}\
                --ad_out {output.H5AD}\
                --feat_col 1 \
                --remove_zero_features
            """
        )

# initialize & cache the **raw** counts as an anndata file for easier loading later
## Removes barcodes for which there are no molecules detected [`--remove_zero_features`]
rule cache_preQC_h5ad_kb:
    input:
        BCS = '{OUTDIR}/{SAMPLE}/kb/{RECIPE}/raw/output.barcodes.txt.gz',
        GENES = '{OUTDIR}/{SAMPLE}/kb/{RECIPE}/raw/output.genes.txt.gz',
        MAT = '{OUTDIR}/{SAMPLE}/kb/{RECIPE}/raw/output.mtx.gz',
        BB_map = lambda wildcards: BB_DICT[wildcards.SAMPLE]
    output:
        H5AD = "{OUTDIR}/{SAMPLE}/kb/{RECIPE}/raw/output.h5ad"
    # params:
    #     var_names = "gene_symbols" # scanpy.read_10x_mtx()
    threads:
        1
    run:
        shell(
            f"""
            python scripts/py/cache_mtx_to_h5ad.py \
                --mat_in {input.MAT} \
                --feat_in {input.GENES} \
                --bc_in {input.BCS} \
                --bb_map {input.BB_map}\
                --ad_out {output.H5AD}\
                --feat_col 0 \
                --remove_zero_features
            """
        )

rule cache_preQC_h5ad_kbpython:
    input:
        BCS = '{OUTDIR}/{SAMPLE}/kbpython/{RECIPE}/counts_unfiltered/cells_x_genes.barcodes.txt.gz',
        GENES = '{OUTDIR}/{SAMPLE}/kbpython/{RECIPE}/counts_unfiltered/cells_x_genes.genes.txt.gz',
        MAT = '{OUTDIR}/{SAMPLE}/kbpython/{RECIPE}/counts_unfiltered/cells_x_genes.mtx.gz',
        BB_map = lambda wildcards: BB_DICT[wildcards.SAMPLE]
    output:
        H5AD = "{OUTDIR}/{SAMPLE}/kbpython/{RECIPE}/counts_unfiltered/output.h5ad"
    # params:
    #     var_names = "gene_symbols" # scanpy.read_10x_mtx()
    threads:
        1
    run:
        shell(
            f"""
            python scripts/py/cache_mtx_to_h5ad.py \
                --mat_in {input.MAT} \
                --feat_in {input.GENES} \
                --bc_in {input.BCS} \
                --bb_map {input.BB_map}\
                --ad_out {output.H5AD}\
                --feat_col 0 \
                --remove_zero_features
            """
        )

# initialize & cache the **raw** counts as an anndata file for easier loading later
## Removes barcodes for which there are no molecules detected [`--remove_zero_features`]
rule cache_preQC_h5ad_miRNA:
    input:
        BCS = '{OUTDIR}/{SAMPLE}/miRNA/{RECIPE}/raw/barcodes.tsv.gz',
        GENES = '{OUTDIR}/{SAMPLE}/miRNA/{RECIPE}/{RECIPE}/raw/features.tsv.gz',
        MAT = '{OUTDIR}/{SAMPLE}/miRNA/{RECIPE}/raw/matrix.mtx.gz',
        BB_map = lambda wildcards: BB_DICT[wildcards.SAMPLE]
    output:
        H5AD = "{OUTDIR}/{SAMPLE}/miRNA/{RECIPE}/raw/output.h5ad"
    threads:
        1
    run:
        shell(
            f"""
            python scripts/py/cache_mtx_to_h5ad.py \
                --mat_in {input.MAT} \
                --feat_in {input.GENES} \
                --bc_in {input.BCS} \
                --bb_map {input.BB_map}\
                --ad_out {output.H5AD}\
                --feat_col 0 \
                --remove_zero_features
            """
        )



# initialize & cache the **raw** counts as an anndata file for easier loading later
## Removes barcodes for which there are no molecules detected [`--remove_zero_features`]
rule cache_preQC_h5ad_piRNA:
    input:
        BCS = '{OUTDIR}/{SAMPLE}/piRNA/{RECIPE}/raw/barcodes.tsv.gz',
        GENES = '{OUTDIR}/{SAMPLE}/piRNA/{RECIPE}/raw/features.tsv.gz',
        MAT = '{OUTDIR}/{SAMPLE}/piRNA/{RECIPE}/raw/matrix.mtx.gz',
        BB_map = lambda wildcards: BB_DICT[wildcards.SAMPLE]
    output:
        H5AD = "{OUTDIR}/{SAMPLE}/piRNA/{RECIPE}/raw/output.h5ad"
    threads:
        1
    run:
        shell(
            f"""
            python scripts/py/cache_mtx_to_h5ad.py \
                --mat_in {input.MAT} \
                --feat_in {input.GENES} \
                --bc_in {input.BCS} \
                --bb_map {input.BB_map}\
                --ad_out {output.H5AD}\
                --feat_col 0 \
                --remove_zero_features
            """
        )

# initialize & cache the **raw** counts as an anndata file for easier loading later
## Removes barcodes for which there are no molecules detected [`--remove_zero_features`]
rule ont_cache_preQC_h5ad_minimap2:
    input:
        COUNTS = '{OUTDIR}/{SAMPLE}/ont/minimap2/{RECIPE}/raw/umitools_counts.tsv.gz',
        BB_map = lambda wildcards: BB_DICT[wildcards.SAMPLE] #TODO Adjust to match barcode handling schemas
    output:
        H5AD = "{OUTDIR}/{SAMPLE}/ont/minimap2/{RECIPE}/raw/output.h5ad"
    threads:
        1
    run:
        shell(
            f"""
            python scripts/py/cache_umitools_h5ad.py \
                {input.COUNTS} \
                {output.H5AD}
            """
                # --bb_map {input.BB_map}\ #TODO
        )