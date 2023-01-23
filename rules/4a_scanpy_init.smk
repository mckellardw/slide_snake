# initialize & cache the **filtered** counts as an anndata file for easier loading later
rule cache_preQC_h5ad_filtered:
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
rule cache_preQC_h5ad_raw:
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