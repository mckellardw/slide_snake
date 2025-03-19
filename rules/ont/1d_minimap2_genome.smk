# Convert gtf to junction bed for minimap2 alignment
rule ont_1d_genome_generate_junction_bed:
    input:
        GTF=lambda wildcards: SAMPLE_SHEET["genes_gtf"][wildcards.SAMPLE],
    output:
        JUNC_BED="{OUTDIR}/{SAMPLE}/ont/minimap2/{RECIPE}/junctions.bed",
    log:
        err="{OUTDIR}/{SAMPLE}/ont/minimap2/{RECIPE}/gff2bed.err",
    conda:
        f"{workflow.basedir}/envs/minimap2.yml"
    shell:
        """
        mkdir -p $(dirname {output.JUNC_BED})

        paftools.js gff2bed \
            -j {input.GTF} \
        > {output.JUNC_BED} \
        2> {log.err}
        """


# Align w/ minimap2
## minimap2 docs - https://lh3.github.io/minimap2/minimap2.html
rule ont_1d_genome_align_minimap2_genome:
    input:
        FQ=lambda w: get_fqs(w, return_type="list", mode="ONT")[1],
        JUNC_BED="{OUTDIR}/{SAMPLE}/ont/minimap2/{RECIPE}/junctions.bed",
    output:
        SAM_TMP=temp("{OUTDIR}/{SAMPLE}/ont/minimap2/{RECIPE}/tmp.sam"),
    params:
        EXTRA_FLAGS=lambda wildcards: RECIPE_SHEET["mm2_extra"][wildcards.RECIPE],
        REF=lambda wildcards: SAMPLE_SHEET["genome_fa"][wildcards.SAMPLE],
        # JUNC_BED=lambda wildcards: SAMPLE_SHEET["mm2_junc_bed"][wildcards.SAMPLE],
    log:
        log="{OUTDIR}/{SAMPLE}/ont/minimap2/{RECIPE}/minimap2.log",
    resources:
        mem="128G",
    threads: config["CORES"]
    conda:
        f"{workflow.basedir}/envs/minimap2.yml"
    shell:
        """
        mkdir -p $(dirname {output.SAM_TMP})

        echo "Genome reference:   {params.REF}" > {log.log} 
        echo "Extra flags:        {params.EXTRA_FLAGS}" >> {log.log} 
        echo "" >> {log.log} 

        minimap2 -ax splice \
            -uf \
            --MD \
            -t {threads} \
            --junc-bed {input.JUNC_BED} \
            {params.EXTRA_FLAGS} {params.REF} \
            {input.FQ} \
        2>> {log.log} \
        > {output.SAM_TMP}
        """


# Sort and compresss minimap2 output
rule ont_1d_genome_sort_compress_output:
    input:
        SAM="{OUTDIR}/{SAMPLE}/ont/minimap2/{RECIPE}/tmp.sam",
    output:
        BAM=temp("{OUTDIR}/{SAMPLE}/ont/minimap2/{RECIPE}/sorted.bam"),
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


# Assign feature (transcript ID) to each alignment
rule ont_1d_genome_featureCounts:
    input:
        BAM="{OUTDIR}/{SAMPLE}/ont/minimap2/{RECIPE}/sorted.bam",
        BAI="{OUTDIR}/{SAMPLE}/ont/minimap2/{RECIPE}/sorted.bam.bai",
    output:
        TSV="{OUTDIR}/{SAMPLE}/ont/minimap2/{RECIPE}/sorted.bam.featureCounts",
        FEAT="{OUTDIR}/{SAMPLE}/ont/minimap2/{RECIPE}/featureCounts.tsv",
    params:
        GTF=lambda wildcards: SAMPLE_SHEET["genes_gtf"][wildcards.SAMPLE],
        EXTRA_FLAGS=lambda wildcards: RECIPE_SHEET["featureCounts_extra"][
            wildcards.RECIPE
        ],
        MIN_TEMPLATE_LENGTH=10,
        MAX_TEMPLATE_LENGTH=10000,
    log:
        log="{OUTDIR}/{SAMPLE}/ont/minimap2/{RECIPE}/featureCounts.log",
    resources:
        mem="32G",
    threads: 1  # long reads can only run single-threaded
    conda:
        f"{workflow.basedir}/envs/minimap2.yml"
    shell:
        """
        featureCounts \
            -a {params.GTF} \
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
        |& tee {log.log}
        """


# -t 'transcript' \
# -g 'transcript_id' \


# --donotsort \
# −−sortReadsByCoordinates \

# TODO?
# rule salmon_quant:
#     input:
#         BAM="{OUTDIR}/{SAMPLE}/ont/minimap2/{RECIPE}/sorted_cb.bam",
#         BAI="{OUTDIR}/{SAMPLE}/ont/minimap2/{RECIPE}/sorted_cb.bam.bai",
#         INDEX="{PATH_TO_SALMON_INDEX}"
#     output:
#         QUANT="{OUTDIR}/{SAMPLE}/ont/minimap2/{RECIPE}/quant.sf"
#     params:
#         LIBTYPE = "A", # Automatic library type detection
#         VALIDATE_MAPPINGS = True, # Validate mappings
#         GC_BIAS = True, # Correct GC bias
#         NUM_GIBBS_SAMPLES = 20 # Number of Gibbs samples
#     log:
#         log = "{OUTDIR}/{SAMPLE}/ont/minimap2/{RECIPE}/salmon_quant.log"
#     threads:
#         config["CORES"]
#     conda:
#          f"{workflow.basedir}/envs/salmon.yml"
#     shell:
#         """
#         salmon quant -i {input.INDEX} -l {params.LIBTYPE} -p {threads} \
#         --validateMappings {params.VALIDATE_MAPPINGS} --gcBias {params.GC_BIAS} \
#         --numGibbsSamples {params.NUM_GIBBS_SAMPLES} -o {output.QUANT} \
#         -1 {input.BAM} 2> {log.log}
#         """


# Add gene tag (GN) to bam...
rule ont_1d_genome_add_featureCounts_to_bam:
    input:
        BAM="{OUTDIR}/{SAMPLE}/ont/minimap2/{RECIPE}/sorted.bam",
        BAI="{OUTDIR}/{SAMPLE}/ont/minimap2/{RECIPE}/sorted.bam.bai",
        TSV="{OUTDIR}/{SAMPLE}/ont/minimap2/{RECIPE}/sorted.bam.featureCounts",
    output:
        BAM="{OUTDIR}/{SAMPLE}/ont/minimap2/{RECIPE}/sorted_gn.bam",
    params:
        READ_ID_COLUMN=0,
        TAG="GN",  # corrected barcode tag
        TAG_COLUMN=3,
    log:
        log="{OUTDIR}/{SAMPLE}/ont/minimap2/{RECIPE}/tsv2tag_1_GN.log",
        err="{OUTDIR}/{SAMPLE}/ont/minimap2/{RECIPE}/tsv2tag_1_GN.err",
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


# Add CB to gene-tagged .bam
rule ont_1d_genome_add_corrected_barcodes:
    input:
        BAM="{OUTDIR}/{SAMPLE}/ont/minimap2/{RECIPE}/sorted_gn.bam",
        TSV="{OUTDIR}/{SAMPLE}/ont/barcodes_umis/{RECIPE}/barcodes_corrected.tsv",
    output:
        BAM="{OUTDIR}/{SAMPLE}/ont/minimap2/{RECIPE}/sorted_gn_cb.bam",
    params:
        READ_ID_COLUMN=0,
        BARCODE_TAG="CB",
        BARCODE_TSV_COLUMN=1,
    log:
        log="{OUTDIR}/{SAMPLE}/ont/minimap2/{RECIPE}/tsv2tag_2_CB.log",
        err="{OUTDIR}/{SAMPLE}/ont/minimap2/{RECIPE}/tsv2tag_2_CB.err",
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
rule ont_1d_genome_add_umis:
    input:
        BAM="{OUTDIR}/{SAMPLE}/ont/minimap2/{RECIPE}/sorted_gn_cb.bam",
        TSV="{OUTDIR}/{SAMPLE}/ont/barcodes_umis/{RECIPE}/barcodes_filtered.tsv",
    output:
        BAM="{OUTDIR}/{SAMPLE}/ont/minimap2/{RECIPE}/sorted_gn_cb_ub.bam",
    params:
        READ_ID_COLUMN=0,
        UMI_TSV_COLUMN=-1,  # last column
        UMI_TAG="UR",  # uncorrected UMI
    log:
        log="{OUTDIR}/{SAMPLE}/ont/minimap2/{RECIPE}/tsv2tag_3_UR.log",
        err="{OUTDIR}/{SAMPLE}/ont/minimap2/{RECIPE}/tsv2tag_3_UR.err",
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
rule ont_1d_genome_filter_bam_empty_tags:
    input:
        BAM="{OUTDIR}/{SAMPLE}/ont/minimap2/{RECIPE}/sorted_gn_cb_ub.bam",
        # BAI="{OUTDIR}/{SAMPLE}/ont/minimap2/{RECIPE}/sorted_gn_cb_ub.bam.bai",
    output:
        BAM="{OUTDIR}/{SAMPLE}/ont/minimap2/{RECIPE}/sorted_filtered_gn_cb_ub.bam",
    params:
        CELL_TAG="CB",  # uncorrected = CR; corrected = CB
        GENE_TAG="GN",  # GN XS
        UMI_TAG="UR",  # uncorrected = UR; corrected = UB
    resources:
        mem="16G",
    threads: 1
    shell:
        """
        samtools view -h {input.BAM} \
        | awk -v tag={params.CELL_TAG} -f scripts/awk/bam_filterEmptyTag.awk \
        | awk -v tag={params.GENE_TAG} -f scripts/awk/bam_filterEmptyTag.awk \
        | awk -v tag={params.UMI_TAG} -f scripts/awk/bam_filterEmptyTag.awk \
        | samtools view -b \
        > {output.BAM}
        """


# Split BAM by strand
rule ont_1d_genome_split_bam_by_strand:
    input:
        BAM="{OUTDIR}/{SAMPLE}/ont/minimap2/{RECIPE}/sorted_filtered_gn_cb_ub.bam",
    output:
        BAM_POS="{OUTDIR}/{SAMPLE}/ont/minimap2/{RECIPE}/sorted_filtered_gn_cb_ub_pos.bam",
        BAM_NEG="{OUTDIR}/{SAMPLE}/ont/minimap2/{RECIPE}/sorted_filtered_gn_cb_ub_neg.bam",
    resources:
        mem="8G",
    threads: 1
    shell:
        """
        samtools view -bh -F 0x10 {input.BAM} > {output.BAM_POS}
        samtools view -bh -f 0x10 {input.BAM} > {output.BAM_NEG}
        """


# Generate count matrix w/ umi-tools
rule ont_1d_genome_umitools_count:
    input:
        BAM="{OUTDIR}/{SAMPLE}/ont/minimap2/{RECIPE}/sorted_filtered_gn_cb_ub.bam",
        BAI="{OUTDIR}/{SAMPLE}/ont/minimap2/{RECIPE}/sorted_filtered_gn_cb_ub.bam.bai",
    output:
        COUNTS="{OUTDIR}/{SAMPLE}/ont/minimap2/{RECIPE}/raw/umitools_counts.tsv.gz",
    params:
        CELL_TAG="CB",  # uncorrected = CR
        GENE_TAG="GN",  #GN XS
        UMI_TAG="UR",
    log:
        log="{OUTDIR}/{SAMPLE}/ont/misc_logs/{RECIPE}/1d_umitools_count.log",
        err="{OUTDIR}/{SAMPLE}/ont/misc_logs/{RECIPE}/1d_minimap2_umitools_count.err",
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
rule ont_1d_genome_counts_to_sparse:
    input:
        COUNTS="{OUTDIR}/{SAMPLE}/ont/minimap2/{RECIPE}/raw/umitools_counts.tsv.gz",
    output:
        BCS="{OUTDIR}/{SAMPLE}/ont/minimap2/{RECIPE}/raw/barcodes.tsv.gz",
        FEATS="{OUTDIR}/{SAMPLE}/ont/minimap2/{RECIPE}/raw/features.tsv.gz",
        COUNTS="{OUTDIR}/{SAMPLE}/ont/minimap2/{RECIPE}/raw/matrix.mtx.gz",
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
rule ont_1d_genome_cache_h5ad_minimap2:
    input:
        BCS="{OUTDIR}/{SAMPLE}/ont/minimap2/{RECIPE}/raw/barcodes.tsv.gz",
        FEATS="{OUTDIR}/{SAMPLE}/ont/minimap2/{RECIPE}/raw/features.tsv.gz",
        MAT="{OUTDIR}/{SAMPLE}/ont/minimap2/{RECIPE}/raw/matrix.mtx.gz",
        BC_map=lambda w: get_bc_map(w, mode="ONT"),
        # BC_map="{OUTDIR}/{SAMPLE}/bc/map_underscore.txt",
    output:
        H5AD="{OUTDIR}/{SAMPLE}/ont/minimap2/{RECIPE}/raw/output.h5ad",
        QC_PLOTS="{OUTDIR}/{SAMPLE}/ont/minimap2/{RECIPE}/raw/qc_h5ad.png",
    log:
        log="{OUTDIR}/{SAMPLE}/ont/minimap2/{RECIPE}/raw/cache_h5ad.log",
        err="{OUTDIR}/{SAMPLE}/ont/minimap2/{RECIPE}/raw/cache_h5ad.err",
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
            --feat_col 1 0 \
            --remove_zero_features \
            --plot_qc \
            --qc_plot_file {output.QC_PLOTS} \
        1> {log.log} \
        2> {log.err}
        """


# make Seurat object with spatial coordinates
rule ont_1d_genome_cache_seurat_minimap2:
    input:
        BCS="{OUTDIR}/{SAMPLE}/ont/minimap2/{RECIPE}/raw/barcodes.tsv.gz",
        FEATS="{OUTDIR}/{SAMPLE}/ont/minimap2/{RECIPE}/raw/features.tsv.gz",
        MAT="{OUTDIR}/{SAMPLE}/ont/minimap2/{RECIPE}/raw/matrix.mtx.gz",
        BC_map=lambda w: get_bc_map(w, mode="ONT"),
    output:
        SEURAT="{OUTDIR}/{SAMPLE}/ont/minimap2/{RECIPE}/raw/output.rds",
    log:
        log="{OUTDIR}/{SAMPLE}/ont/minimap2/{RECIPE}/raw/cache_seurat.log",
        err="{OUTDIR}/{SAMPLE}/ont/minimap2/{RECIPE}/raw/cache_seurat.err",
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
