# Align w/ ultra
## ultra docs - https://github.com/ksahlin/ultra?tab=readme-ov-file#USAGE

# rule ont_2d_ultra_index_genome:
#     input:
#         FA=lambda wildcards: SAMPLE_SHEET["genome_fa"][wildcards.SAMPLE],
#     output:
#         SAM_TMP=temp("{OUTDIR}/{SAMPLE}/ont/ultra/{RECIPE}/tmp.sam"),
#     params:
#         FA=lambda wildcards: SAMPLE_SHEET["genome_fa"][wildcards.SAMPLE],
#         JUNC_BED=lambda wildcards: SAMPLE_SHEET["mm2_junc_bed"][wildcards.SAMPLE],
#     log:
#         log="{OUTDIR}/REFERENCES/ont/ultra/{RECIPE}/logs/ultra.log",
#         err="{OUTDIR}/REFERENCES/ont/ultra/{RECIPE}/logs/ultra.err",
#     resources:
#         mem="128G",
#     threads: config["CORES"]
#     conda:
#         f"{workflow.basedir}/envs/ultra.yml"
#     shell:
#         """
#         uLTRA index \
#             {input.FA} \
#             {input.GTF} \
#             $(dirname {output.FILE}) \
#             --disable_infer \
#         1> {log.log} \
#         2> {log.err}
#         """


# Sort input gtf
## Prevents this bug- https://github.com/ksahlin/ultra/issues/24
rule ont_1f_sort_gtf:
    input:
        GTF=lambda wildcards: SAMPLE_SHEET["genes_gtf"][wildcards.SAMPLE],
    output:
        GTF=temp("{OUTDIR}/{SAMPLE}/ont/ultra/{RECIPE}/sorted.gtf"),
    resources:
        mem="128G",
    threads: 1
    conda:
        f"{workflow.basedir}/envs/ultra.yml"
    shell:
        """
        gff3sort.pl {input.GTF} > {output.GTF}
        """


rule ont_2d_ultra_pipeline_genome:
    input:
        FQ=lambda w: get_fqs(w, return_type="list", mode="ONT")[1],
        GTF="{OUTDIR}/{SAMPLE}/ont/ultra/{RECIPE}/sorted.gtf",
    output:
        SAM=temp("{OUTDIR}/{SAMPLE}/ont/ultra/{RECIPE}/reads.sam"),
        DB="{OUTDIR}/{SAMPLE}/ont/ultra/{RECIPE}/database.db",
    params:
        EXTRA_FLAGS=lambda wildcards: RECIPE_SHEET["ultra_extra"][wildcards.RECIPE],
        FA=lambda wildcards: SAMPLE_SHEET["genome_fa"][wildcards.SAMPLE],
        # GTF=lambda wildcards: SAMPLE_SHEET["genes_gtf"][wildcards.SAMPLE],
    log:
        log="{OUTDIR}/{SAMPLE}/ont/ultra/{RECIPE}/logs/ultra.log",
        err="{OUTDIR}/{SAMPLE}/ont/ultra/{RECIPE}/logs/ultra.err",
    resources:
        mem="128G",
    # threads: config["CORES"]
    conda:
        f"{workflow.basedir}/envs/ultra.yml"
    shell:
        """
        uLTRA pipeline \
            {params.FA} \
            {input.GTF} \
            {input.FQ} \
            $(dirname {output.SAM}) \
            --ont \
            --t {threads} {params.EXTRA_FLAGS} \
            --disable_infer \
        1> {log.log} \
        2> {log.err}
        """


# Sort and compresss ultra output
rule ont_2d_ultra_sort_compress_output:
    input:
        SAM="{OUTDIR}/{SAMPLE}/ont/ultra/{RECIPE}/reads.sam",
    output:
        BAM=temp("{OUTDIR}/{SAMPLE}/ont/ultra/{RECIPE}/sorted.bam"),
    params:
        REF=lambda wildcards: SAMPLE_SHEET["genome_fa"][wildcards.SAMPLE],
    resources:
        mem="16G",
    threads: 1
    shell:
        """
        samtools sort --reference {params.REF} \
            -O BAM \
            -o {output.BAM} \
            {input.SAM}             
        """


# Add CB to BAM
rule ont_2d_ultra_add_corrected_barcodes:
    input:
        BAM="{OUTDIR}/{SAMPLE}/ont/ultra/{RECIPE}/sorted.bam",
        TSV="{OUTDIR}/{SAMPLE}/ont/barcodes_umis/{RECIPE}/barcodes_corrected.tsv",
    output:
        BAM="{OUTDIR}/{SAMPLE}/ont/ultra/{RECIPE}/sorted_cb.bam",
    params:
        READ_ID_COLUMN=0,
        BARCODE_TAG="CB",  # corrected barcode
        BARCODE_TSV_COLUMN=1,
    log:
        log="{OUTDIR}/{SAMPLE}/ont/ultra/{RECIPE}/logs/tsv2tag_1_CB.log",
        err="{OUTDIR}/{SAMPLE}/ont/ultra/{RECIPE}/logs/tsv2tag_1_CB.err",
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
            --tagColumns {params.BARCODE_TSV_COLUMN} \
            --tags {params.BARCODE_TAG} \
        1> {log.log} \
        2> {log.err}
        """


# Add UMI (UR) to barcoded BAM
rule ont_2d_ultra_add_umis:
    input:
        BAM="{OUTDIR}/{SAMPLE}/ont/ultra/{RECIPE}/sorted_cb.bam",
        TSV="{OUTDIR}/{SAMPLE}/ont/barcodes_umis/{RECIPE}/barcodes_filtered.tsv",
    output:
        BAM="{OUTDIR}/{SAMPLE}/ont/ultra/{RECIPE}/sorted_cb_ub.bam",
    params:
        READ_ID_COLUMN=0,
        UMI_TSV_COLUMN=-1,  # last column
        UMI_TAG="UR",  # uncorrected UMI
    log:
        log="{OUTDIR}/{SAMPLE}/ont/ultra/{RECIPE}/logs/tsv2tag_2_UR.log",
        err="{OUTDIR}/{SAMPLE}/ont/ultra/{RECIPE}/logs/tsv2tag_2_UR.err",
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
            --tagColumns {params.UMI_TSV_COLUMN} \
            --tags {params.UMI_TAG} \
        1> {log.log} \
        2> {log.err}
        """


# Generate count matrix w/ umi-tools
rule ont_2d_ultra_filter_bam_empty_tags:
    input:
        BAM="{OUTDIR}/{SAMPLE}/ont/ultra/{RECIPE}/sorted_cb_ub.bam",
        # BAI="{OUTDIR}/{SAMPLE}/ont/ultra/{RECIPE}/sorted_cb_ub.bam.bai",
    output:
        BAM="{OUTDIR}/{SAMPLE}/ont/ultra/{RECIPE}/sorted_filtered_cb_ub.bam",
    params:
        CELL_TAG="CB",  # uncorrected = CR; corrected = CB
        # GENE_TAG="GN",  # GN XS
        UMI_TAG="UR",  # uncorrected = UR; corrected = UB
    resources:
        mem="16G",
    threads: 1
    shell:
        """
        samtools view -h {input.BAM} \
        | awk -v tag={params.CELL_TAG} -f scripts/awk/bam_filterEmptyTag.awk \
        | awk -v tag={params.UMI_TAG} -f scripts/awk/bam_filterEmptyTag.awk \
        | samtools view -b \
        > {output.BAM}
        """


# Assign feature (transcript ID) to each alignment
rule ont_2d_ultra_featureCounts:
    input:
        BAM="{OUTDIR}/{SAMPLE}/ont/ultra/{RECIPE}/sorted_filtered_cb_ub.bam",
        BAI="{OUTDIR}/{SAMPLE}/ont/ultra/{RECIPE}/sorted_filtered_cb_ub.bam.bai",
        GTF="{OUTDIR}/{SAMPLE}/ont/ultra/{RECIPE}/sorted.gtf",
    output:
        TSV="{OUTDIR}/{SAMPLE}/ont/ultra/{RECIPE}/featureCounts/sorted_filtered_cb_ub.bam.featureCounts",
        FEAT="{OUTDIR}/{SAMPLE}/ont/ultra/{RECIPE}/featureCounts/featureCounts.tsv",
    params:
        GTF=lambda wildcards: SAMPLE_SHEET["genes_gtf"][wildcards.SAMPLE],
        EXTRA_FLAGS=lambda wildcards: RECIPE_SHEET["featureCounts_extra"][
            wildcards.RECIPE
        ],
        MIN_TEMPLATE_LENGTH=10,
        MAX_TEMPLATE_LENGTH=10000,
    log:
        log="{OUTDIR}/{SAMPLE}/ont/ultra/{RECIPE}/featureCounts/featureCounts.log",
        err="{OUTDIR}/{SAMPLE}/ont/ultra/{RECIPE}/featureCounts/featureCounts.err",
    resources:
        mem="32G",
    threads: 1  # long reads can only run single-threaded
    conda:
        f"{workflow.basedir}/envs/minimap2.yml"
    shell:
        """
        mkdir -p $(dirname {output.TSV})
        featureCounts \
            -a {input.GTF} \
            -o {output.FEAT} \
            -L \
            -s 1 \
            -f \
            -d {params.MIN_TEMPLATE_LENGTH} \
            -D {params.MAX_TEMPLATE_LENGTH} \
            -t 'gene' \
            -g 'gene_id' \
            -T {threads} \
            -R CORE {params.EXTRA_FLAGS} \
            {input.BAM} \
        1> {log.log} \
        2> {log.err} \
        """


# Add gene tag (GN) to BAM
rule ont_2d_ultra_add_featureCounts_to_bam:
    input:
        BAM="{OUTDIR}/{SAMPLE}/ont/ultra/{RECIPE}/sorted_filtered_cb_ub.bam",
        BAI="{OUTDIR}/{SAMPLE}/ont/ultra/{RECIPE}/sorted_filtered_cb_ub.bam.bai",
        TSV="{OUTDIR}/{SAMPLE}/ont/ultra/{RECIPE}/featureCounts/sorted_filtered_cb_ub.bam.featureCounts",
    output:
        BAM="{OUTDIR}/{SAMPLE}/ont/ultra/{RECIPE}/sorted_filtered_cb_ub_gn.bam",
    params:
        READ_ID_COLUMN=0,
        TAG="GN",  # gene tag
        TAG_COLUMN=3,
    log:
        log="{OUTDIR}/{SAMPLE}/ont/ultra/{RECIPE}/logs/tsv2tag_3_GN.log",
        err="{OUTDIR}/{SAMPLE}/ont/ultra/{RECIPE}/logs/tsv2tag_3_GN.err",
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
rule ont_2d_ultra_umitools_count:
    input:
        BAM="{OUTDIR}/{SAMPLE}/ont/ultra/{RECIPE}/sorted_filtered_cb_ub_gn.bam",
        BAI="{OUTDIR}/{SAMPLE}/ont/ultra/{RECIPE}/sorted_filtered_cb_ub_gn.bam.bai",
    output:
        COUNTS="{OUTDIR}/{SAMPLE}/ont/ultra/{RECIPE}/raw/umitools_counts.tsv.gz",
    params:
        CELL_TAG="CB",  # uncorrected = CR
        GENE_TAG="GN",  #GN XS
        UMI_TAG="UR",
    log:
        log="{OUTDIR}/{SAMPLE}/ont/ultra/{RECIPE}/logs/umitools_count.log",
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
            -S {output.COUNTS} 
        """


# Convert long-format counts from umi_tools to market-matrix format (.mtx)
rule ont_2d_ultra_counts_to_sparse:
    input:
        COUNTS="{OUTDIR}/{SAMPLE}/ont/ultra/{RECIPE}/raw/umitools_counts.tsv.gz",
    output:
        BCS="{OUTDIR}/{SAMPLE}/ont/ultra/{RECIPE}/raw/barcodes.tsv.gz",
        FEATS="{OUTDIR}/{SAMPLE}/ont/ultra/{RECIPE}/raw/features.tsv.gz",
        COUNTS="{OUTDIR}/{SAMPLE}/ont/ultra/{RECIPE}/raw/matrix.mtx.gz",
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
rule ont_2d_ultra_cache_h5ad:
    input:
        BCS="{OUTDIR}/{SAMPLE}/ont/ultra/{RECIPE}/raw/barcodes.tsv.gz",
        FEATS="{OUTDIR}/{SAMPLE}/ont/ultra/{RECIPE}/raw/features.tsv.gz",
        MAT="{OUTDIR}/{SAMPLE}/ont/ultra/{RECIPE}/raw/matrix.mtx.gz",
        BC_map=lambda w: get_bc_map(w, mode="ONT"),
        # BC_map="{OUTDIR}/{SAMPLE}/bc/map_underscore.txt",
    output:
        H5AD="{OUTDIR}/{SAMPLE}/ont/ultra/{RECIPE}/raw/output.h5ad",
        QC_PLOTS="{OUTDIR}/{SAMPLE}/ont/ultra/{RECIPE}/raw/qc_plots.png",
    params:
        GTF=lambda w: SAMPLE_SHEET["genes_gtf"][w.SAMPLE],
        GTF_FEATURE_TYPE="gene",  # feature type in gtf to use 
        GTF_ID="gene_id",  # gtf attribute used to match var_names in adata
    log:
        log="{OUTDIR}/{SAMPLE}/ont/ultra/{RECIPE}/logs/cache_h5ad.log",
        err="{OUTDIR}/{SAMPLE}/ont/ultra/{RECIPE}/logs/cache_h5ad.err",
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
            --qc_plot_file {output.QC_PLOTS} \
            --gtf_file {params.GTF} \
            --gtf_feature_type {params.GTF_FEATURE_TYPE} \
            --gtf_id {params.GTF_ID} \
        1> {log.log} \
        2> {log.err}
        """


# make Seurat object with spatial coordinates
rule ont_2d_ultra_cache_seurat:
    input:
        BCS="{OUTDIR}/{SAMPLE}/ont/ultra/{RECIPE}/raw/barcodes.tsv.gz",
        FEATS="{OUTDIR}/{SAMPLE}/ont/ultra/{RECIPE}/raw/features.tsv.gz",
        MAT="{OUTDIR}/{SAMPLE}/ont/ultra/{RECIPE}/raw/matrix.mtx.gz",
        BC_map=lambda w: get_bc_map(w, mode="ONT"),
    output:
        SEURAT="{OUTDIR}/{SAMPLE}/ont/ultra/{RECIPE}/raw/output.rds",
    params:
        FEAT_COL=1,  # 1=gene_name, 2=gene_id
        TRANSPOSE="True",
    log:
        log="{OUTDIR}/{SAMPLE}/ont/ultra/{RECIPE}/logs/cache_seurat.log",
        err="{OUTDIR}/{SAMPLE}/ont/ultra/{RECIPE}/logs/cache_seurat.err",
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
