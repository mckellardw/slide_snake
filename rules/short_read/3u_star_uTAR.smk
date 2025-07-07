# ported from https://github.com/mckellardw/uTAR_lite
# Link to paper: https://www.nature.com/articles/s41467-021-22496-3


# Filter BAM to remove reads without GN tag
## Deduped BAM file already has reads with no CB/UB tags removed
rule ilmn_3u_filter_noGN:
    input:
        BAM="{OUTDIR}/{SAMPLE}/short_read/STARsolo/{RECIPE}/Aligned.sortedByCoord.out.dedup.bam",
    output:
        BAM="{OUTDIR}/{SAMPLE}/short_read/STARsolo/{RECIPE}/Aligned.sortedByCoord.noGN.dedup.bam",
    params:
        TAG="GN",
    shell:
        """
        samtools view -h {input.BAM} \
        | awk -v tag={params.TAG} -f scripts/awk/bam_keepEmptyTag.awk \
        | samtools view -b -o {output.BAM}
        """


# Run HMM
rule ilmn_3u_calcHMMbed:
    input:
        BAM="{OUTDIR}/{SAMPLE}/short_read/STARsolo/{RECIPE}/Aligned.sortedByCoord.noGN.dedup.bam",
        BAI="{OUTDIR}/{SAMPLE}/short_read/STARsolo/{RECIPE}/Aligned.sortedByCoord.noGN.dedup.bam.bai",
    output:
        BED=temp("{OUTDIR}/{SAMPLE}/short_read/STARsolo/{RECIPE}/TAR/TAR_raw.bed.gz"),
    threads: config["CORES"]
    resources:
        mem_mb=65536,  # 64GB in MB
    params:
        # MEM="64G",
        MERGEBP=100,  # default 100 (Note- window size across genome is 50bp)
        THRESH=10000000,  # default 10000000 (Note- number of reads to use for HMM)
    log:
        log="{OUTDIR}/{SAMPLE}/short_read/STARsolo/{RECIPE}/TAR/run_hmm.log",
        err="{OUTDIR}/{SAMPLE}/short_read/STARsolo/{RECIPE}/TAR/run_hmm.err",
    conda:
        f"{workflow.basedir}/envs/utar.yml"
    shell:
        """
        bash scripts/bash/bam_uTAR_HMM.sh \
            --bam {input.BAM} \
            --threads {threads} \
            --mem {resources.mem_mb} \
            --mergebp {params.MERGEBP} \
            --thresh {params.THRESH} \
            --outdir $(dirname {output.BED}) \
        1> {log.log} \
        2> {log.err}
        """


rule ilmn_3u_filter_out_aTARs:
    input:
        BED="{OUTDIR}/{SAMPLE}/short_read/STARsolo/{RECIPE}/TAR/TAR_raw.bed.gz",
    output:
        BED="{OUTDIR}/{SAMPLE}/short_read/STARsolo/{RECIPE}/TAR/uTAR.bed.gz",
    params:
        GTF=lambda wildcards: SAMPLE_SHEET["genes_gtf"][wildcards.SAMPLE],
    resources:
        mem_mb=8000,
    conda:
        f"{workflow.basedir}/envs/ucsc.yml"
    shell:
        """
        gtfToGenePred {params.GTF} stdout \
        | genePredToBed stdin stdout \
        | bedtools intersect -v -a {input.BED} -b stdin \
        | bgzip > {output.BED}
        """


# Convert BED to GTF
rule ilmn_3u_bed_to_gtf:
    input:
        BEDGZ="{OUTDIR}/{SAMPLE}/short_read/STARsolo/{RECIPE}/TAR/uTAR.bed.gz",
    output:
        GTF="{OUTDIR}/{SAMPLE}/short_read/STARsolo/{RECIPE}/TAR/uTAR.gtf",
    conda:
        f"{workflow.basedir}/envs/ucsc.yml"
    shell:
        """
        zcat {input.BEDGZ} \
        | bedToGenePred stdin stdout \
        | genePredToGtf file stdin {output.GTF}
        """


# Label .bam file with each HMM feature
rule ilmn_3u_tagReads:
    input:
        BAM="{OUTDIR}/{SAMPLE}/short_read/STARsolo/{RECIPE}/Aligned.sortedByCoord.noGN.dedup.bam",
        TAR_GTF="{OUTDIR}/{SAMPLE}/short_read/STARsolo/{RECIPE}/TAR/uTAR.gtf",
    output:
        DIR=directory(
            "{OUTDIR}/{SAMPLE}/short_read/STARsolo/{RECIPE}/TAR/gene_assigned"
        ),
        BAM=temp(
            "{OUTDIR}/{SAMPLE}/short_read/STARsolo/{RECIPE}/TAR/Aligned.sortedByCoord.noGN.dedup.bam.featureCounts.bam"
        ),
    threads: config["CORES"]
    resources:
        mem_mb=16000,
    log:
        log="{OUTDIR}/{SAMPLE}/short_read/STARsolo/{RECIPE}/TAR/featureCounts.log",
        err="{OUTDIR}/{SAMPLE}/short_read/STARsolo/{RECIPE}/TAR/featureCounts.err",
    conda:
        f"{workflow.basedir}/envs/umi_tools.yml"
    shell:
        """
        featureCounts \
            -T {threads} \
            -t exon \
            -g gene_id \
            -a {input.TAR_GTF} \
            --largestOverlap \
            --readExtension5 0 \
            --readExtension3 0 \
            -s 1 \
            -M \
            -o {output.DIR} \
            -R BAM \
            {input.BAM} \
        1> {log.log} \
        2> {log.err}
        """


rule ilmn_3u_sort_index_tagged_bam:
    input:
        BAM="{OUTDIR}/{SAMPLE}/short_read/STARsolo/{RECIPE}/TAR/Aligned.sortedByCoord.noGN.dedup.bam.featureCounts.bam",
    output:
        BAM="{OUTDIR}/{SAMPLE}/short_read/STARsolo/{RECIPE}/TAR/TAR_tagged.bam",
    threads: config["CORES"]
    shell:
        """
        samtools sort -@ {threads} {input.BAM} -o {output.BAM}
        """


# Get counts matrix for HMM-annotated features
rule ilmn_3u_extract_HMM_expression:
    input:
        BAM="{OUTDIR}/{SAMPLE}/short_read/STARsolo/{RECIPE}/TAR/TAR_tagged.bam",
        BAI="{OUTDIR}/{SAMPLE}/short_read/STARsolo/{RECIPE}/TAR/TAR_tagged.bam.bai",
    output:
        COUNT_MTX="{OUTDIR}/{SAMPLE}/short_read/STARsolo/{RECIPE}/TAR/counts.tsv.gz",
    params:
        CELL_TAG="CB",  # uncorrected = CR
        GENE_TAG="XT",  #GN XS
        UMI_TAG="UB",
        STATUS_TAG="XS",
    log:
        log="{OUTDIR}/{SAMPLE}/short_read/STARsolo/{RECIPE}/TAR/umitools_count.log",
        err="{OUTDIR}/{SAMPLE}/short_read/STARsolo/{RECIPE}/TAR/umitools_count.err",
    conda:
        f"{workflow.basedir}/envs/umi_tools.yml"
    shell:
        """
        umi_tools count --extract-umi-method=tag \
            --per-cell \
            --per-gene \
            --assigned-status-tag={params.CELL_TAG} \
            --cell-tag={params.CELL_TAG} \
            --gene-tag={params.GENE_TAG}  \
            --umi-tag={params.UMI_TAG}  \
            --log={log.log} \
            --stdin {input.BAM} \
            --stdout {output.COUNT_MTX} \
        2> {log.err}
        """
        # --multimapping-detection-method=NH \


# Convert the long-format matrix (from umi_tools) to an .mtx file
rule ilmn_3u_counts_long2mtx:
    input:
        TSV="{OUTDIR}/{SAMPLE}/short_read/STARsolo/{RECIPE}/TAR/counts.tsv.gz",
    output:
        MTX="{OUTDIR}/{SAMPLE}/short_read/STARsolo/{RECIPE}/TAR/uTAR.mtx",
        GENES="{OUTDIR}/{SAMPLE}/short_read/STARsolo/{RECIPE}/TAR/uTAR_genes.tsv.gz",
        CELLS="{OUTDIR}/{SAMPLE}/short_read/STARsolo/{RECIPE}/TAR/uTAR_cells.tsv.gz",
    conda:
        f"{workflow.basedir}/envs/scanpy.yml"
    shell:
        """
        python scripts/py/long2mtx.py \
            --umitools_tsv {input.TSV} \
            --out_mat {output.MTX} \
            --output-format mtx
        """


# Compress the count matrix
rule ilmn_3u_gzip_counts:
    input:
        COUNT_MTX="{OUTDIR}/{SAMPLE}/short_read/STARsolo/{RECIPE}/TAR/uTAR.mtx",
    output:
        COUNT_MTX="{OUTDIR}/{SAMPLE}/short_read/STARsolo/{RECIPE}/TAR/uTAR.mtx.gz",
    threads: config["CORES"]
    shell:
        """
        pigz -p{threads} {input.COUNT_MTX}
        """


# Generate QC plots for the final matrix
rule ilmn_3u_plot_qc:
    input:
        MTX="{OUTDIR}/{SAMPLE}/short_read/STARsolo/{RECIPE}/TAR/uTAR.mtx.gz",
        GENES="{OUTDIR}/{SAMPLE}/short_read/STARsolo/{RECIPE}/TAR/uTAR_genes.tsv.gz",
        CELLS="{OUTDIR}/{SAMPLE}/short_read/STARsolo/{RECIPE}/TAR/uTAR_cells.tsv.gz",
    output:
        PNG="{OUTDIR}/{SAMPLE}/short_read/STARsolo/{RECIPE}/TAR/qc_plots.png",
    conda:
        f"{workflow.basedir}/envs/seurat.yml"
    shell:
        """
        Rscript scripts/R/matrix_qc_summary.R \
            --matrix {input.MTX} \
            --genes {input.GENES} \
            --cells {input.CELLS} \
            --output {output.PNG}
        """
