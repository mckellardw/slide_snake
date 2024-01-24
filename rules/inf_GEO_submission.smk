# Prep fastq files for GEO upload
rule prep_GEO_fqs:
    input:
        BCS = '{OUTDIR}/{sample}/miRNA/raw/barcodes.tsv.gz',
        GENES = '{OUTDIR}/{sample}/miRNA/raw/features.tsv.gz',
        MAT = '{OUTDIR}/{sample}/miRNA/raw/matrix.mtx.gz',
        BB_map = lambda wildcards: BB_DICT[wildcards.sample]
    output:
        H5AD = "{OUTDIR}/{sample}/miRNA/raw/output.h5ad"
    threads:
        1
    run:
        shell(
            f"""
            echo "TODO"
            """
        )

# Prep count matrices for GEO
rule prep_GEO_counts:
    input:
        BCS = '{OUTDIR}/{sample}/miRNA/raw/barcodes.tsv.gz',
        GENES = '{OUTDIR}/{sample}/miRNA/raw/features.tsv.gz',
        MAT = '{OUTDIR}/{sample}/miRNA/raw/matrix.mtx.gz',
        BB_map = lambda wildcards: BB_DICT[wildcards.sample]
    output:
        H5AD = "{OUTDIR}/{sample}/miRNA/raw/output.h5ad"
    threads:
        1
    run:
        shell(
            f"""
            echo "TODO"
            """
        )