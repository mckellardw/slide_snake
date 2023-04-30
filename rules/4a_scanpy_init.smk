# initialize & cache the **filtered** counts as an anndata file for easier loading later
rule cache_preQC_h5ad_STAR_filtered:
    input:
        GENEFULLMAT = "{DATADIR}/align_out/{sample}/STARsolo/Solo.out/GeneFull/filtered/matrix.mtx.gz"
    output:
        H5AD = "{DATADIR}/align_out/{sample}/STARsolo/Solo.out/GeneFull/filtered/matrix.h5ad"
    params:
        var_names = "gene_symbols" # scanpy.read_10x_mtx()
    threads:
        1
    run:
        shell(
            f"""
            python scripts/cache_h5ad.py {DATADIR}/align_out/{wildcards.sample}/STARsolo/Solo.out/GeneFull/filtered {output.H5AD} {params.var_names}
            """
        )

# initialize & cache the **raw** counts as an anndata file for easier loading later
rule cache_preQC_h5ad_STAR_raw:
    input:
        GENEFULLMAT = "{DATADIR}/align_out/{sample}/STARsolo/Solo.out/GeneFull/raw/UniqueAndMult-EM.mtx.gz"
    output:
        H5AD = "{DATADIR}/align_out/{sample}/STARsolo/Solo.out/GeneFull/raw/UniqueAndMultEM.h5ad"
    params:
        var_names = "gene_symbols" # scanpy.read_10x_mtx()
    threads:
        1
    run:
        shell(
            f"""
            python scripts/cache_h5ad.py {DATADIR}/align_out/{wildcards.sample}/STARsolo/Solo.out/GeneFull/raw {output.H5AD} {params.var_names}
            """
        )

# initialize & cache the **raw** counts as an anndata file for easier loading later
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
            {input.MAT} \
            {input.GENES} \
            {input.BCS} \
            {input.BB_map} \
            {output.H5AD}
            """
        )