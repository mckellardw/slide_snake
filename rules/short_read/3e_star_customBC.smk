# Add custom barcode calling results to STAR-aligned .bam and generate count matrices


# Use STAR-barcoded .bam file & filter:
#   - Don't filter by length, as the longest hairpin miRNA (in mouse) is 147nt
#   - Remove reads missing a cell barcode ("CB" tag) or UMI ("UB" tag)
#   - Remove all aligment info, but keep tags
#   - Remove STAR alignment score ("AS" tag)
#   - Clip the 3' end of the read, to remove non-templated additions
# rule ilmn_3N_STARcustom_clearCBTag:
#     input:
#         BAM="{OUTDIR}/{SAMPLE}/STARsolo/Aligned.sortedByCoord.out.bam",
#     output:
#         BAM="{OUTDIR}/{SAMPLE}/STARcustom/Aligned.sortedByCoord.out.bam",
#     params:
#         UNCORRECTED_BARCODE_TAG="UB",
#         CORRECTED_BARCODE_TAG="CB",
#     threads: 1
#     shell:
#         """
#         mkdir -p $(dirname {output.BAM})

#         samtools view {input.BAM} \
#         | awk -v tag={params.UNCORRECTED_BARCODE_TAG} -f scripts/awk/bam_filterEmptyTag.awk - \
#         | awk -v tag={params.CORRECTED_BARCODE_TAG} -f scripts/awk/bam_filterEmptyTag.awk - \
#         | samtools view -bS \
#         > {output.BAM}
#         """


#############################################
## STARcustom: Align reads using STAR (2-Pass; R2 only)
#############################################
rule ilmn_3e_STARcustom_firstPass:
    input:
        FQS=lambda w: get_fqs(w, return_type="list", mode="ILMN"),
    output:
        SJ="{OUTDIR}/{SAMPLE}/short_read/STARcustom/{RECIPE}/firstPass/_STARpass1/SJ.out.tab",
    params:
        STAR_REF=lambda w: get_STAR_ref(w),
        STAR_EXTRA=lambda w: get_STAR_extra_params(w)["STAR_extra"],
    log:
        log="{OUTDIR}/{SAMPLE}/short_read/STARcustom/{RECIPE}/firstPass/pass1.log",
        err="{OUTDIR}/{SAMPLE}/short_read/STARcustom/{RECIPE}/firstPass/pass1.err",
    resources:
        mem=megabytes2bytes(config["MEMLIMIT_MB"]),
        time="2:00:00",
    threads: config["CORES"]
    conda:
        f"{workflow.basedir}/envs/star.yml"
    shell:
        """
        mkdir -p $(dirname {output.SJ})
        STAR --runThreadN {threads} \
            --outFileNamePrefix $(dirname {output.SJ})/ \
            --genomeDir {params.STAR_REF} \
            --readFilesCommand zcat \
            --readFilesIn {input.FQS[1]} \
            --outSAMtype BAM Unsorted {params.STAR_EXTRA} \
            --twopassMode Basic \
        1> {log.log} 2> {log.err}
        """


rule ilmn_3e_STARcustom_secondPass:
    input:
        FQS=lambda w: get_fqs(w, return_type="list", mode="ILMN"),
        SJ="{OUTDIR}/{SAMPLE}/short_read/STARcustom/{RECIPE}/firstPass/_STARpass1/SJ.out.tab",
    output:
        BAM="{OUTDIR}/{SAMPLE}/short_read/STARcustom/{RECIPE}/Aligned.sortedByCoord.out.bam",
        UNMAPPED=[
            f"{{OUTDIR}}/{{SAMPLE}}/short_read/STARcustom/{{RECIPE}}/Unmapped.out.mate{READ}"
            for READ in [1, 2]
        ],
    params:
        STAR_REF=lambda w: get_STAR_ref(w),
        outBAMsortingBinsN=100,
        outBAMsortingThreadN=4,
    log:
        log="{OUTDIR}/{SAMPLE}/short_read/STARcustom/{RECIPE}/pass2.log",
        err="{OUTDIR}/{SAMPLE}/short_read/STARcustom/{RECIPE}/pass2.err",
    resources:
        mem=config["MEMLIMIT"],
        time="2:00:00",
    threads: config["CORES"]
    conda:
        f"{workflow.basedir}/envs/star.yml"
    shell:
        """
        mkdir -p $(dirname {output.BAM})
        STAR --readFilesCommand zcat \
            --readFilesIn {input.FQS[1]} \
            --genomeDir {params.STAR_REF} \
            --runThreadN {threads} \
            --limitBAMsortRAM={resources.mem} \
            --outBAMsortingBinsN {params.outBAMsortingBinsN} \
            --outBAMsortingThreadN {params.outBAMsortingThreadN} \
            --outFileNamePrefix $(dirname {output.BAM})/ \
            --outSAMtype BAM SortedByCoordinate \
            --outSAMattributes NH HI nM AS CR UR CB UB GX GN sS sQ sM \
            --outReadsUnmapped Fastx --outSAMunmapped Within KeepPairs \
            --sjdbFileChrStartEnd {input.SJ} \
        1> {log.log} \
        2> {log.err}
        """


rule ilmn_3e_gzip_unmapped_fastq:
    input:
        FQ_R1="{OUTDIR}/{SAMPLE}/short_read/STARcustom/{RECIPE}/Unmapped.out.mate1",
        FQ_R2="{OUTDIR}/{SAMPLE}/short_read/STARcustom/{RECIPE}/Unmapped.out.mate2",
    output:
        FQ_R1="{OUTDIR}/{SAMPLE}/short_read/STARcustom/{RECIPE}/Unmapped.out.mate1.fq.gz",
        FQ_R2="{OUTDIR}/{SAMPLE}/short_read/STARcustom/{RECIPE}/Unmapped.out.mate2.fq.gz",
    threads: config["CORES"]
    shell:
        """
        pigz -p {threads} -c {input.FQ_R1} > {output.FQ_R1}
        pigz -p {threads} -c {input.FQ_R2} > {output.FQ_R2}
        """


# Add CB to gene-tagged .bam
rule ilmn_3N_STARcustom_add_corrected_barcodes:
    input:
        BAM="{OUTDIR}/{SAMPLE}/short_read/STARcustom/{RECIPE}/Aligned.sortedByCoord.out.bam",
        TSV="{OUTDIR}/{SAMPLE}/short_read/barcodes_umis/{RECIPE}/read_barcodes_corrected.tsv",
    output:
        BAM="{OUTDIR}/{SAMPLE}/short_read/STARcustom/{RECIPE}/sorted_gn_cb.bam",
    params:
        READ_ID_COLUMN=0,
        BARCODE_TAG="CB",
        BARCODE_TSV_COLUMN=1,
    log:
        log="{OUTDIR}/{SAMPLE}/short_read/STARcustom/{RECIPE}/tsv2tag_2_CB.log",
        err="{OUTDIR}/{SAMPLE}/short_read/STARcustom/{RECIPE}/tsv2tag_2_CB.err",
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


# Split BAM by strand
rule ilmn_3N_STARcustom_split_bam_by_strand:
    input:
        BAM="{OUTDIR}/{SAMPLE}/short_read/STARcustom/{RECIPE}/sorted_filtered_gn_cb_ub.bam",
    output:
        BAM_POS="{OUTDIR}/{SAMPLE}/short_read/STARcustom/{RECIPE}/sorted_filtered_gn_cb_ub_pos.bam",
        BAM_NEG="{OUTDIR}/{SAMPLE}/short_read/STARcustom/{RECIPE}/sorted_filtered_gn_cb_ub_neg.bam",
    resources:
        mem="8G",
    threads: 1
    shell:
        """
        samtools view -bh -F 0x10 {input.BAM} > {output.BAM_POS}
        samtools view -bh -f 0x10 {input.BAM} > {output.BAM_NEG}
        """


# Generate count matrix w/ umi-tools
rule ilmn_3N_STARcustom_umitools_count:
    input:
        BAM="{OUTDIR}/{SAMPLE}/short_read/STARcustom/{RECIPE}/sorted_filtered_gn_cb_ub.bam",
        BAI="{OUTDIR}/{SAMPLE}/short_read/STARcustom/{RECIPE}/sorted_filtered_gn_cb_ub.bam.bai",
    output:
        COUNTS="{OUTDIR}/{SAMPLE}/short_read/STARcustom/{RECIPE}/raw/umitools_counts.tsv.gz",
    params:
        CELL_TAG="CB",  # uncorrected = CR
        GENE_TAG="GN",  #GN XS
        UMI_TAG="UR",
    log:
        log="{OUTDIR}/{SAMPLE}/short_read/misc_logs/{RECIPE}/1d_umitools_count.log",
        err="{OUTDIR}/{SAMPLE}/short_read/misc_logs/{RECIPE}/1d_minimap2_umitools_count.err",
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
rule ilmn_3N_STARcustom_counts_to_sparse:
    input:
        COUNTS="{OUTDIR}/{SAMPLE}/short_read/STARcustom/{RECIPE}/raw/umitools_counts.tsv.gz",
    output:
        BCS="{OUTDIR}/{SAMPLE}/short_read/STARcustom/{RECIPE}/raw/barcodes.tsv.gz",
        FEATS="{OUTDIR}/{SAMPLE}/short_read/STARcustom/{RECIPE}/raw/features.tsv.gz",
        COUNTS="{OUTDIR}/{SAMPLE}/short_read/STARcustom/{RECIPE}/raw/matrix.mtx.gz",
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
rule ilmn_3N_STARcustom_cache_preQC_h5ad_minimap2:
    input:
        BCS="{OUTDIR}/{SAMPLE}/short_read/STARcustom/{RECIPE}/raw/barcodes.tsv.gz",
        FEATS="{OUTDIR}/{SAMPLE}/short_read/STARcustom/{RECIPE}/raw/features.tsv.gz",
        MAT="{OUTDIR}/{SAMPLE}/short_read/STARcustom/{RECIPE}/raw/matrix.mtx.gz",
        BC_map=lambda w: get_bc_map(w, mode="short_read"),
    output:
        H5AD="{OUTDIR}/{SAMPLE}/short_read/STARcustom/{RECIPE}/raw/output.h5ad",
        QC_PLOTS="{OUTDIR}/{SAMPLE}/short_read/STARcustom/{RECIPE}/raw/qc_plots.png",
    log:
        log="{OUTDIR}/{SAMPLE}/short_read/misc_logs/{RECIPE}/1d_minimap2_cache.log",
        err="{OUTDIR}/{SAMPLE}/short_read/misc_logs/{RECIPE}/1d_minimap2_cache.err",
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
        1> {log.log} \
        2> {log.err}
        """
