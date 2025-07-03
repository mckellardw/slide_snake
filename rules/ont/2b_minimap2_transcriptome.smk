# Align w/ minimap2
## minimap2 docs - https://lh3.github.io/minimap2/minimap2.html
## see oarfish requirements - https://github.com/COMBINE-lab/oarfish?tab=readme-ov-file#choosing-minimap2-alignment-options
rule ont_2b_txome_align_minimap2_transcriptome:
    input:
        FQ=lambda w: get_fqs(w, return_type="list", mode="ONT")[1],
    output:
        BAM=temp("{OUTDIR}/{SAMPLE}/ont/minimap2_txome/{RECIPE}/aligned.bam"),
    params:
        EXTRA_FLAGS=lambda wildcards: RECIPE_SHEET["mm2_extra"][wildcards.RECIPE],
        REF=lambda wildcards: SAMPLE_SHEET["cdna_fa"][wildcards.SAMPLE],
    log:
        log="{OUTDIR}/{SAMPLE}/ont/minimap2_txome/{RECIPE}/minimap2.log",
    resources:
        mem="128G",
    threads: config["CORES"]
    conda:
        f"{workflow.basedir}/envs/minimap2.yml"
    shell:
        """
        mkdir -p $(dirname {output.BAM})

        echo "Genome reference:   {params.REF}" > {log.log} 
        echo "Junction reference: {params.REF}" >> {log.log} 
        echo "Extra flags:        {params.EXTRA_FLAGS}" >> {log.log} 
        echo "" >> {log.log} 

        minimap2 \
            -ax map-ont \
            --eqx \
            -N 100 \
            -t {threads} \
            {params.REF} \
            {input.FQ} \
        | samtools view -b \
        > {output.BAM} \
        2>> {log.log}
        """


# Add CB to gene-tagged .bam
rule ont_2b_txome_add_corrected_barcodes:
    input:
        BAM="{OUTDIR}/{SAMPLE}/ont/minimap2_txome/{RECIPE}/aligned.bam",
        TSV="{OUTDIR}/{SAMPLE}/ont/barcodes_umis/{RECIPE}/barcodes_corrected.tsv",
    output:
        BAM="{OUTDIR}/{SAMPLE}/ont/minimap2_txome/{RECIPE}/aligned_cb.bam",
    params:
        READ_ID_COLUMN=0,
        BARCODE_TAG="CB",  # corrected barcode
        BARCODE_TSV_COLUMN=1,
    log:
        log="{OUTDIR}/{SAMPLE}/ont/minimap2_txome/{RECIPE}/tsv2tag_1_CB.log",
        err="{OUTDIR}/{SAMPLE}/ont/minimap2_txome/{RECIPE}/tsv2tag_1_CB.err",
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


# Add UMI (UR) to barcoded & gene-tagged .bam
rule ont_2b_txome_add_umis:
    input:
        BAM="{OUTDIR}/{SAMPLE}/ont/minimap2_txome/{RECIPE}/aligned_cb.bam",
        TSV="{OUTDIR}/{SAMPLE}/ont/barcodes_umis/{RECIPE}/barcodes_filtered.tsv",
    output:
        BAM="{OUTDIR}/{SAMPLE}/ont/minimap2_txome/{RECIPE}/aligned_cb_ub.bam",
    params:
        READ_ID_COLUMN=0,
        UMI_TSV_COLUMN=-1,  # last column
        UMI_TAG="UR",  # uncorrected UMI
    log:
        log="{OUTDIR}/{SAMPLE}/ont/minimap2_txome/{RECIPE}/tsv2tag_2_UR.log",
        err="{OUTDIR}/{SAMPLE}/ont/minimap2_txome/{RECIPE}/tsv2tag_2_UR.err",
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
rule ont_2b_txome_filter_bam_empty_tags:
    input:
        BAM="{OUTDIR}/{SAMPLE}/ont/minimap2_txome/{RECIPE}/aligned_cb_ub.bam",
        # BAI="{OUTDIR}/{SAMPLE}/ont/minimap2_txome/{RECIPE}/aligned_gn_cb_ub.bam.bai",
    output:
        BAM="{OUTDIR}/{SAMPLE}/ont/minimap2_txome/{RECIPE}/aligned_filtered_cb_ub_xb.bam",
    params:
        CELL_TAG="CB",  # uncorrected = CR; corrected = CB
        GENE_TAG="GN",  # GN XS
        UMI_TAG="UR",  # uncorrected = UR; corrected = UB
        COMBINED_TAGS="XB",  # combined tag
    resources:
        mem="16G",
    threads: 1
    shell:
        """
        samtools view -h {input.BAM} \
        | awk -v tag={params.CELL_TAG} -f scripts/awk/bam_filterEmptyTag.awk \
        | awk -v tag={params.UMI_TAG} -f scripts/awk/bam_filterEmptyTag.awk \
        | awk -v tags={params.CELL_TAG},{params.UMI_TAG} -v outtag={params.COMBINED_TAGS} -f bam_combineTags.awk \
        | samtools view -b \
        > {output.BAM}
        """


# Sort aligned_filtered_cb_ub.bam by CB tag
## see more oarfish requirements here - https://github.com/COMBINE-lab/oarfish?tab=readme-ov-file#notes-about-single-cell-mode
rule ont_2b_txome_sort_by_cb:
    input:
        BAM="{OUTDIR}/{SAMPLE}/ont/minimap2_txome/{RECIPE}/aligned_filtered_cb_ub_xb.bam",
    output:
        BAM="{OUTDIR}/{SAMPLE}/ont/minimap2_txome/{RECIPE}/aligned_filtered_sorted_cb_ub_xb.bam",
    params:
        BC_TAG="CB",  # corrected barcode tag
    log:
        err="{OUTDIR}/{SAMPLE}/ont/minimap2_txome/{RECIPE}/sort_by_cb.err",
    resources:
        mem="16G",
    threads: 1
    conda:
        f"{workflow.basedir}/envs/minimap2.yml"
    shell:
        """
        samtools sort -t {params.BC_TAG} -o {output.BAM} {input.BAM} 2> {log.err}
        """

# Deduplicate BAM file based on XB tag using umi_tools
rule ont_2b_txome_dedup_by_xb:
    input:
        BAM="{OUTDIR}/{SAMPLE}/ont/minimap2_txome/{RECIPE}/aligned_filtered_sorted_cb_ub_xb.bam",
    output:
        BAM="{OUTDIR}/{SAMPLE}/ont/minimap2_txome/{RECIPE}/aligned_filtered_sorted_cb_ub_xb_dedup.bam",
    log:
        log="{OUTDIR}/{SAMPLE}/ont/minimap2_txome/{RECIPE}/dedup_by_xb.log",
        err="{OUTDIR}/{SAMPLE}/ont/minimap2_txome/{RECIPE}/dedup_by_xb.err",
    params:
        TAG="XB",
    threads: 1
    resources:
        mem="16G",
    conda:
        f"{workflow.basedir}/envs/umi_tools.yml"
    shell:
        """
        umi_tools dedup \
            --extract-umi-method=tag \
            --umi-tag={params.TAG} \
            -I {input.BAM} \
            -S {output.BAM} \
            1> {log.log} \
            2> {log.err}
        """

# Run oarfish alignment mode transcript quantification
# github: https://github.com/COMBINE-lab/oarfish
# TODO- --short-quant?
rule ont_2b_txome_oarfish_quant:
    input:
        BAM="{OUTDIR}/{SAMPLE}/ont/minimap2_txome/{RECIPE}/aligned_filtered_sorted_cb_ub_xb.bam",
        # BAI="{OUTDIR}/{SAMPLE}/ont/minimap2_txome/{RECIPE}/aligned_filtered_sorted_cb_ub.bam.bai",
    output:
        DIR=directory("{OUTDIR}/{SAMPLE}/ont/minimap2_txome/{RECIPE}/oarfish"),
        JSON="{OUTDIR}/{SAMPLE}/ont/minimap2_txome/{RECIPE}/oarfish/P.meta_info.json",
        QUANT="{OUTDIR}/{SAMPLE}/ont/minimap2_txome/{RECIPE}/oarfish/P.quant",
        PQ="{OUTDIR}/{SAMPLE}/ont/minimap2_txome/{RECIPE}/oarfish/P.infreps.pq",
        AMBIG="{OUTDIR}/{SAMPLE}/ont/minimap2_txome/{RECIPE}/oarfish/P.ambig_info.tsv",
    log:
        log="{OUTDIR}/{SAMPLE}/ont/minimap2_txome/{RECIPE}/oarfish/oarfish.log",
        err="{OUTDIR}/{SAMPLE}/ont/minimap2_txome/{RECIPE}/oarfish/oarfish.err",
    resources:
        mem="32G",
    threads: 16
    conda:
        f"{workflow.basedir}/envs/oarfish.yml"
    shell:
        """
        mkdir -p {output.DIR}

        oarfish \
            --alignments {input.BAM} \
            --output {output.DIR}/P \
            --threads {threads} \
            --filter-group no-filters \
            --model-coverage \
            --single-cell \
            --verbose \
        1> >(awk -f scripts/awk/remove_ansi.awk > {log.log}) \
        2> >(awk -f scripts/awk/remove_ansi.awk > {log.err}
        """


# Code snippet to remove the annoying ANSI codes from the log file:
## sed -r 's/\x1B\[[0-9;]*[mK]//g' oarfish.err > oarfish_cleaned.err


# Generate count matrix w/ umi-tools
rule ont_2b_txome_umitools_count:
    input:
        BAM="{OUTDIR}/{SAMPLE}/ont/minimap2_txome/{RECIPE}/aligned_filtered_gn_cb_ub.bam",
        BAI="{OUTDIR}/{SAMPLE}/ont/minimap2_txome/{RECIPE}/aligned_filtered_gn_cb_ub.bam.bai",
    output:
        COUNTS="{OUTDIR}/{SAMPLE}/ont/minimap2_txome/{RECIPE}/raw/umitools_counts.tsv.gz",
    params:
        CELL_TAG="CB",  # uncorrected = CR
        GENE_TAG="GN",  #GN XS
        UMI_TAG="UR",
    log:
        log="{OUTDIR}/{SAMPLE}/ont/minimap2_txome/{RECIPE}/umitools_count.log",
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
rule ont_2b_txome_counts_to_sparse:
    input:
        COUNTS="{OUTDIR}/{SAMPLE}/ont/minimap2_txome/{RECIPE}/raw/umitools_counts.tsv.gz",
    output:
        BCS="{OUTDIR}/{SAMPLE}/ont/minimap2_txome/{RECIPE}/raw/barcodes.tsv.gz",
        FEATS="{OUTDIR}/{SAMPLE}/ont/minimap2_txome/{RECIPE}/raw/features.tsv.gz",
        COUNTS="{OUTDIR}/{SAMPLE}/ont/minimap2_txome/{RECIPE}/raw/matrix.mtx.gz",
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
rule ont_2b_txome_cache_preQC_h5ad_minimap2:
    input:
        BCS="{OUTDIR}/{SAMPLE}/ont/minimap2_txome/{RECIPE}/raw/barcodes.tsv.gz",
        FEATS="{OUTDIR}/{SAMPLE}/ont/minimap2_txome/{RECIPE}/raw/features.tsv.gz",
        MAT="{OUTDIR}/{SAMPLE}/ont/minimap2_txome/{RECIPE}/raw/matrix.mtx.gz",
        BC_map=lambda w: get_bc_map(w, mode="ONT"),
        # BC_map="{OUTDIR}/{SAMPLE}/bc/map_underscore.txt",
    output:
        H5AD="{OUTDIR}/{SAMPLE}/ont/minimap2_txome/{RECIPE}/raw/output.h5ad",
        QC_PLOTS="{OUTDIR}/{SAMPLE}/ont/minimap2_txome/{RECIPE}/raw/qc_plots.png",
    log:
        log="{OUTDIR}/{SAMPLE}/ont/minimap2_txome/{RECIPE}/raw/cache_h5ad.log",
        err="{OUTDIR}/{SAMPLE}/ont/minimap2_txome/{RECIPE}/raw/cache_h5ad.err",
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
        1> {log.log}
        """


# make Seurat object with spatial coordinates
rule ont_2b_txome_cache_preQC_seurat_minimap2:
    input:
        BCS="{OUTDIR}/{SAMPLE}/ont/minimap2_txome/{RECIPE}/raw/barcodes.tsv.gz",
        FEATS="{OUTDIR}/{SAMPLE}/ont/minimap2_txome/{RECIPE}/raw/features.tsv.gz",
        MAT="{OUTDIR}/{SAMPLE}/ont/minimap2_txome/{RECIPE}/raw/matrix.mtx.gz",
        BC_map=lambda w: get_bc_map(w, mode="ONT"),
    output:
        SEURAT="{OUTDIR}/{SAMPLE}/ont/minimap2_txome/{RECIPE}/raw/output.rds",
    log:
        log="{OUTDIR}/{SAMPLE}/ont/minimap2_txome/{RECIPE}/raw/cache.log",
        err="{OUTDIR}/{SAMPLE}/ont/minimap2_txome/{RECIPE}/raw/cache.err",
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
