# isoquant documentation - https://ablab.github.io/IsoQuant/
# TODO- dedup bam before isoquant
rule ont_2e_isoquant:
    input:
        BAM="{OUTDIR}/{SAMPLE}/ont/{ALIGNER}/{RECIPE}/sorted_filtered_cb_ub_gn.bam",
        BAI="{OUTDIR}/{SAMPLE}/ont/{ALIGNER}/{RECIPE}/sorted_filtered_cb_ub_gn.bam.bai",
    output:
        MAT_TX="{OUTDIR}/{SAMPLE}/ont/{ALIGNER}/{RECIPE}/isoquant/isoquant.transcript_grouped_counts.tsv",
        MAT_GN="{OUTDIR}/{SAMPLE}/ont/{ALIGNER}/{RECIPE}/isoquant/isoquant.gene_grouped_counts.tsv",
        GTF="{OUTDIR}/{SAMPLE}/ont/{ALIGNER}/{RECIPE}/isoquant/isoquant.extended_annotation.gtf",
        TSV="{OUTDIR}/{SAMPLE}/ont/{ALIGNER}/{RECIPE}/isoquant/isoquant.read_assignments.tsv.gz",
    params:
        OUTDIR="{OUTDIR}/{SAMPLE}/ont/{ALIGNER}/{RECIPE}/",
        GENES_GTF=lambda wildcards: SAMPLE_SHEET["genes_gtf"][wildcards.SAMPLE],
        REF=lambda wildcards: SAMPLE_SHEET["genome_fa"][wildcards.SAMPLE],
        CELL_TAG="CB",  # uncorrected = CR; corrected = CB
        # GENE_TAG="GN",  # GN XS
        # UMI_TAG="UR",  # uncorrected = UR; corrected = UB
    log:
        log="{OUTDIR}/{SAMPLE}/ont/{ALIGNER}/{RECIPE}/logs/isoquant.log",
        err="{OUTDIR}/{SAMPLE}/ont/{ALIGNER}/{RECIPE}/logs/isoquant.err",
    threads: config["CORES"]
    conda:
        f"{workflow.basedir}/envs/isoquant.yml"
    shell:
        """
        isoquant.py \
            --reference {params.REF} \
            --genedb {params.GENES_GTF} \
            --bam {input.BAM} \
            --data_type nanopore \
            --stranded forward \
            --polya_requirement never \
            --report_canonical all \
            --check_canonical \
            --gene_quantification all \
            --transcript_quantification all \
            --report_novel_unspliced true \
            --splice_correction_strategy default_ont \
            --complete_genedb \
            --sqanti_output \
            --read_group tag:{params.CELL_TAG} \
            --counts_format matrix \
            --force \
            --threads {threads} \
            --prefix isoquant \
            --output {params.OUTDIR} \
        1> {log.log} \
        2> {log.err}
        """


# Add isoquant gene tag (IG) to bam...
rule ont_2e_add_isoquant_genes_to_bam:
    input:
        BAM="{OUTDIR}/{SAMPLE}/ont/{ALIGNER}/{RECIPE}/sorted_filtered_cb_ub_gn.bam",
        BAI="{OUTDIR}/{SAMPLE}/ont/{ALIGNER}/{RECIPE}/sorted_filtered_cb_ub_gn.bam.bai",
        TSV="{OUTDIR}/{SAMPLE}/ont/{ALIGNER}/{RECIPE}/isoquant/isoquant.read_assignments.tsv.gz",
    output:
        BAM="{OUTDIR}/{SAMPLE}/ont/{ALIGNER}/{RECIPE}/sorted_filtered_cb_ub_gn_ig.bam",
    params:
        READ_ID_COLUMN=0,
        TAG="IG",  # corrected barcode tag
        TAG_COLUMN=4,  # column in tsv with gene tag
    log:
        log="{OUTDIR}/{SAMPLE}/ont/{ALIGNER}/{RECIPE}/logs/tsv2tag_4_IG.log",
        err="{OUTDIR}/{SAMPLE}/ont/{ALIGNER}/{RECIPE}/logs/tsv2tag_4_IG.err",
    conda:
        f"{workflow.basedir}/envs/parasail.yml"
    resources:
        mem="16G",
    threads: 1
    shell:
        """
        python scripts/py/tsv2tag.py --in_bam {input.BAM} \
            --in_tsv {input.TSV} \
            --out_bam {output.BAM} \
            --readIDColumn {params.READ_ID_COLUMN} \
            --tagColumns {params.TAG_COLUMN} \
            --tags {params.TAG} \
        1> {log.log} \
        2> {log.err}
        """


# Add isoquant transcript tag (IT) to bam...
rule ont_2e_add_isoquant_transcripts_to_bam:
    input:
        BAM="{OUTDIR}/{SAMPLE}/ont/{ALIGNER}/{RECIPE}/sorted_filtered_cb_ub_gn_ig.bam",
        BAI="{OUTDIR}/{SAMPLE}/ont/{ALIGNER}/{RECIPE}/sorted_filtered_cb_ub_gn_ig.bam.bai",
        TSV="{OUTDIR}/{SAMPLE}/ont/{ALIGNER}/{RECIPE}/isoquant/isoquant.read_assignments.tsv.gz",
    output:
        BAM="{OUTDIR}/{SAMPLE}/ont/{ALIGNER}/{RECIPE}/sorted_filtered_cb_ub_gn_ig_it.bam",
    params:
        READ_ID_COLUMN=0,
        TAG="IT",  # corrected barcode tag
        TAG_COLUMN=3,  # column in tsv with gene tag
    log:
        log="{OUTDIR}/{SAMPLE}/ont/{ALIGNER}/{RECIPE}/logs/tsv2tag_5_IT.log",
        err="{OUTDIR}/{SAMPLE}/ont/{ALIGNER}/{RECIPE}/logs/tsv2tag_5_IT.err",
    conda:
        f"{workflow.basedir}/envs/parasail.yml"
    resources:
        mem="16G",
    threads: 1
    shell:
        """
        python scripts/py/tsv2tag.py --in_bam {input.BAM} \
            --in_tsv {input.TSV} \
            --out_bam {output.BAM} \
            --readIDColumn {params.READ_ID_COLUMN} \
            --tagColumns {params.TAG_COLUMN} \
            --tags {params.TAG} \
        1> {log.log} \
        2> {log.err}
        """


# Generate count matrix w/ umi-tools
rule ont_2e_umitools_count:
    input:
        BAM="{OUTDIR}/{SAMPLE}/ont/{ALIGNER}/{RECIPE}/sorted_filtered_cb_ub_gn_ig_it.bam",
        BAI="{OUTDIR}/{SAMPLE}/ont/{ALIGNER}/{RECIPE}/sorted_filtered_cb_ub_gn_ig_it.bam.bai",
    output:
        COUNTS="{OUTDIR}/{SAMPLE}/ont/{ALIGNER}/{RECIPE}/isoquant/umitools_counts.tsv.gz",
    params:
        CELL_TAG="CB",  # uncorrected = CR; corrected = CB
        GENE_TAG="IT",  # isoquant transcript tag
        # GENE_TAG="IG",  # isoquant gene tag
        UMI_TAG="UR",  # uncorrected = UR; corrected = UB
    log:
        log="{OUTDIR}/{SAMPLE}/ont/{ALIGNER}/{RECIPE}/logs/umitools_isoquant.log",
        err="{OUTDIR}/{SAMPLE}/ont/{ALIGNER}/{RECIPE}/logs/umitools_isoquant.err",
    resources:
        mem="16G",
    threads: 1
    conda:
        f"{workflow.basedir}/envs/umi_tools.yml"
    shell:
        """
        umi_tools count --extract-umi-method=tag \
            --per-gene \
            --per-cell \
            --cell-tag={params.CELL_TAG} \
            --gene-tag={params.GENE_TAG}  \
            --umi-tag={params.UMI_TAG}  \
            --log={log.log} \
            -I {input.BAM} \
            -S {output.COUNTS} \
        2> {log.err}
        """


# Convert long-format counts from umi_tools to market-matrix format (.mtx)
rule ont_2e_counts_to_sparse:
    input:
        COUNTS="{OUTDIR}/{SAMPLE}/ont/{ALIGNER}/{RECIPE}/isoquant/umitools_counts.tsv.gz",
    output:
        BCS="{OUTDIR}/{SAMPLE}/ont/{ALIGNER}/{RECIPE}/isoquant/barcodes.tsv.gz",
        FEATS="{OUTDIR}/{SAMPLE}/ont/{ALIGNER}/{RECIPE}/isoquant/features.tsv.gz",
        COUNTS="{OUTDIR}/{SAMPLE}/ont/{ALIGNER}/{RECIPE}/isoquant/matrix.mtx.gz",
    resources:
        mem="16G",
    threads: 1
    conda:
        f"{workflow.basedir}/envs/scanpy.yml"
    shell:
        """
        mkdir -p $(dirname {output.COUNTS})
        python scripts/py/long2mtx.py {input.COUNTS} $(dirname {output.COUNTS})
        """


# make anndata object with spatial coordinates
rule ont_2e_cache_h5ad:
    input:
        BCS="{OUTDIR}/{SAMPLE}/ont/{ALIGNER}/{RECIPE}/isoquant/barcodes.tsv.gz",
        FEATS="{OUTDIR}/{SAMPLE}/ont/{ALIGNER}/{RECIPE}/isoquant/features.tsv.gz",
        MAT="{OUTDIR}/{SAMPLE}/ont/{ALIGNER}/{RECIPE}/isoquant/matrix.mtx.gz",
        BC_map=lambda w: get_bc_map(w, mode="ONT"),
    output:
        H5AD="{OUTDIR}/{SAMPLE}/ont/{ALIGNER}/{RECIPE}/isoquant/output.h5ad",
        QC_PLOTS="{OUTDIR}/{SAMPLE}/ont/{ALIGNER}/{RECIPE}/isoquant/qc_h5ad.png",
    params:
        GTF=lambda w: SAMPLE_SHEET["genes_gtf"][w.SAMPLE],
        GTF_FEATURE_TYPE="gene",  # feature type in gtf to use 
        GTF_ID="gene_id",  # gtf attribute used to match var_names in adata
    log:
        log="{OUTDIR}/{SAMPLE}/ont/{ALIGNER}/{RECIPE}/logs/cache_isoquant_h5ad.log",
        err="{OUTDIR}/{SAMPLE}/ont/{ALIGNER}/{RECIPE}/logs/cache_isoquant_h5ad.err",
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
            --feat_col 0 \
            --remove_zero_features \
            --plot_qc \
            --qc_plot_file {output.QC_PLOTS} \
            --gtf_file {params.GTF} \
            --gtf_feature_type {params.GTF_FEATURE_TYPE} \
            --gtf_id {params.GTF_ID} \
        1> {log.log} \
        2> {log.err}
        """


# make Seurat object with spatial coordinates
rule ont_2e_cache_seurat:
    input:
        BCS="{OUTDIR}/{SAMPLE}/ont/{ALIGNER}/{RECIPE}/isoquant/barcodes.tsv.gz",
        FEATS="{OUTDIR}/{SAMPLE}/ont/{ALIGNER}/{RECIPE}/isoquant/features.tsv.gz",
        MAT="{OUTDIR}/{SAMPLE}/ont/{ALIGNER}/{RECIPE}/isoquant/matrix.mtx.gz",
        BC_map=lambda w: get_bc_map(w, mode="ONT"),
    output:
        SEURAT="{OUTDIR}/{SAMPLE}/ont/{ALIGNER}/{RECIPE}/isoquant/output.rds",
    params:
        FEAT_COL=1,  # column in features file to use as feature names (R is 1-indexed)
        TRANSPOSE="False",  # whether to transpose the matrix (default is False)
    log:
        log="{OUTDIR}/{SAMPLE}/ont/{ALIGNER}/{RECIPE}/logs/cache_isoquant_seurat.log",
        err="{OUTDIR}/{SAMPLE}/ont/{ALIGNER}/{RECIPE}/logs/cache_isoquant_seurat.err",
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
            --feat_col {params.FEAT_COL} \
            --transpose {params.TRANSPOSE} \
        1> {log.log} \
        2> {log.err}
        """
