# Align w/ minimap2
## minimap2 docs - https://lh3.github.io/minimap2/minimap2.html
rule ont_1d_align_minimap2_genome:
    input:
        # FQ="{OUTDIR}/{SAMPLE}/ont/umitools/{RECIPE}/umi_R2.fq.gz",
        FQ=lambda w: get_fqs(w, return_type="list", mode="ONT")[1],
    output:
        SAM_TMP=temp("{OUTDIR}/{SAMPLE}/ont/minimap2/{RECIPE}/tmp.sam"),
    params:
        EXTRA_FLAGS=lambda wildcards: RECIPE_SHEET["mm2_extra"][wildcards.RECIPE],
        # ref=config["REF_GENOME_FASTA"],
        # chrom_sizes=config["REF_CHROM_SIZES"],
        # bed=config["REF_GENES_BED"],
        # flags=config["RESOURCES_MM2_FLAGS"],
        REF=lambda wildcards: SAMPLE_SHEET["mm2_fa"][wildcards.SAMPLE],
        JUNC_BED=lambda wildcards: SAMPLE_SHEET["mm2_junc_bed"][wildcards.SAMPLE],
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
        echo "Junction reference: {params.REF}" >> {log.log} 
        echo "Extra flags:        {params.EXTRA_FLAGS}" >> {log.log} 
        echo "" >> {log.log} 

        minimap2 -ax splice \
            -uf \
            --MD \
            -t {threads} \
            --junc-bed {params.JUNC_BED} \
            {params.EXTRA_FLAGS} {params.REF} \
            {input.FQ} \
        2>> {log.log} \
        > {output.SAM_TMP}
        """


rule ont_1d_sort_compress_output:
    input:
        SAM="{OUTDIR}/{SAMPLE}/ont/minimap2/{RECIPE}/tmp.sam",
    output:
        # BAM_UNSORT_TMP=temp("{OUTDIR}/{SAMPLE}/ont/tmp_unsort.sam"),
        BAM=temp("{OUTDIR}/{SAMPLE}/ont/minimap2/{RECIPE}/sorted.bam"),
    params:
        REF=lambda wildcards: SAMPLE_SHEET["mm2_fa"][wildcards.SAMPLE],
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


rule ont_1d_featureCounts:
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
            -t 'transcript' \
            -g 'transcript_id' \
            -T {threads} \
            -R CORE {params.EXTRA_FLAGS} \
            {input.BAM} \
        |& tee {log.log}
        """


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
rule ont_1d_add_featureCounts_to_bam:
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
        |& tee {log.log}
        """


# Add CB to gene-tagged .bam
rule ont_1d_add_corrected_barcodes:
    input:
        BAM="{OUTDIR}/{SAMPLE}/ont/minimap2/{RECIPE}/sorted_gn.bam",
        TSV="{OUTDIR}/{SAMPLE}/ont/barcodes_umis/{RECIPE}/read_barcodes_corrected.tsv",
    output:
        BAM="{OUTDIR}/{SAMPLE}/ont/minimap2/{RECIPE}/sorted_gn_cb.bam",
    params:
        READ_ID_COLUMN=0,
        BARCODE_TAG="CB",  # corrected barcode
        BARCODE_TSV_COLUMN=1,
    log:
        log="{OUTDIR}/{SAMPLE}/ont/minimap2/{RECIPE}/tsv2tag_2_CB.log",
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
        |& tee {log.log}
        """


# Add UMI (UR) to barcoded & gene-tagged .bam
rule ont_1d_add_umis:
    input:
        BAM="{OUTDIR}/{SAMPLE}/ont/minimap2/{RECIPE}/sorted_gn_cb.bam",
        TSV="{OUTDIR}/{SAMPLE}/ont/barcodes_umis/{RECIPE}/read_barcodes_filtered.tsv",
    output:
        BAM="{OUTDIR}/{SAMPLE}/ont/minimap2/{RECIPE}/sorted_gn_cb_ub.bam",
    params:
        READ_ID_COLUMN=0,
        UMI_TSV_COLUMN=-1,  # last column
        UMI_TAG="UR",  # uncorrected UMI
    log:
        log="{OUTDIR}/{SAMPLE}/ont/minimap2/{RECIPE}/tsv2tag_3_UR.log",
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
        |& tee {log.log}
        """


# Generate count matrix w/ umi-tools
rule ont_1d_filter_bam_empty_tags:
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


# Generate count matrix w/ umi-tools
rule ont_1d_umitools_count:
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
        log="{OUTDIR}/{SAMPLE}/ont/minimap2/{RECIPE}/umitools_count.log",
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


rule ont_1d_counts_to_sparse:
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
