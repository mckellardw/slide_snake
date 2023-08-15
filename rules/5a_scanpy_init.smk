# initialize & cache the **raw** counts as an anndata file for easier loading later
## Removes barcodes for which there are no molecules detected [`--remove_zero_features`]
rule cache_preQC_h5ad_STAR_raw:
    input:
        BCS = '{OUTDIR}/{sample}/STARsolo/Solo.out/GeneFull/raw/barcodes.tsv.gz',
        GENES = '{OUTDIR}/{sample}/STARsolo/Solo.out/GeneFull/raw/features.tsv.gz',
        MAT = '{OUTDIR}/{sample}/STARsolo/Solo.out/GeneFull/raw/UniqueAndMult-EM.mtx.gz',
        BB_map = lambda wildcards: BB_DICT[wildcards.sample]
    output:
        H5AD = '{OUTDIR}/{sample}/STARsolo/Solo.out/GeneFull/raw/UniqueAndMultEM.h5ad'
    params:
        var_names = "gene_symbols" # scanpy.read_10x_mtx()
    threads:
        1
    run:
        shell(
            f"""
            python scripts/cache_mtx_to_h5ad.py \
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
rule cache_preQC_h5ad_kb_raw:
    input:
        BCS = '{OUTDIR}/{sample}/kb/counts_unfiltered/output.barcodes.txt.gz',
        GENES = '{OUTDIR}/{sample}/kb/counts_unfiltered/output.genes.txt.gz',
        MAT = '{OUTDIR}/{sample}/kb/counts_unfiltered/output.mtx.gz',
        BB_map = lambda wildcards: BB_DICT[wildcards.sample]
    output:
        H5AD = "{OUTDIR}/{sample}/kb/counts_unfiltered/output.h5ad"
    # params:
    #     var_names = "gene_symbols" # scanpy.read_10x_mtx()
    threads:
        1
    run:
        shell(
            f"""
            python scripts/cache_mtx_to_h5ad.py \
            --mat_in {input.MAT} \
            --feat_in {input.GENES} \
            --bc_in {input.BCS} \
            --bb_map {input.BB_map}\
            --ad_out {output.H5AD}\
            --feat_col 0 \
            --remove_zero_features
            """
        )