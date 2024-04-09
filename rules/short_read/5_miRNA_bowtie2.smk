####################################################
# Quantification of (mature) miRNAs

# Summary of Workflow:

# 1) bowtie2_prep_bam: Filter and preprocess the STAR-aligned BAM file to remove
#       reads longer than miRNAs, and clear alignment information while keeping tags.

# 2) bowtie2_miRNA: Align the preprocessed BAM file to the miRNA reference using
#       bowtie2, and save the aligned BAM file.

# 3) umitools_sortAlignedBAM_miRNA: Sort the aligned BAM file using samtools.

# 4) umitools_dedupSortedBAM_miRNA: Deduplicate the sorted BAM file using umi-tools.

# 5) umitools_indexSortedDedupBAM_miRNA: Index the deduplicated BAM file using samtools.

# 6) bowtie2_miRNA_counts: Tag the deduplicated BAM file with chromosome/miRNA information
#       and generate a count matrix using umi-tools.


# bowtie2 documentation: https://bowtie-bio.sourceforge.net/bowtie2/manual.shtml
# pirbase download link - http://bigdata.ibp.ac.cn/piRBase/download.php
#       *Note* - using gold standard miRNA refs ()
####################################################


# Use STAR-barcoded .bam file & filter:
#   - Don't filter by length, as the longest hairpin miRNA (in mouse) is 147nt
#   - Remove reads missing a cell barcode ("CB" tag) or UMI ("UB" tag)
#   - Remove all aligment info, but keep tags
#   - Remove STAR alignment score ("AS" tag)
#   - Clip the 3' end of the read, to remove non-templated additions
rule bowtie2_prep_bam_miRNA:
    input:
        BAM="{OUTDIR}/{SAMPLE}/STARsolo/Aligned.sortedByCoord.dedup.out.bam",
    output:
        BAM=temp("{OUTDIR}/{SAMPLE}/miRNA/{RECIPE}/tmp.bam"),
    params:
        MEMLIMIT=config["MEMLIMIT"],
        MAX_SHORT_READ_LENGTH=config["MAX_SHORT_READ_LENGTH"],
        TRIM_N=4,
        TRIM_OPTION="end",  # option: start, end, both
    threads: config["CORES"]
    # 1
    run:
        shell(
            f"""
            mkdir -p $(dirname {output.BAM})

            {EXEC['SAMTOOLS']} view {input.BAM} \
            | awk -v tag=CB -f scripts/awk/bam_filterEmptyTag.awk - \
            | awk -v tag=UB -f scripts/awk/bam_filterEmptyTag.awk - \
            | awk -f scripts/awk/bam_clearAlignment.awk - \
            | awk -v tag=AS -f scripts/awk/bam_clearTag.awk - \
            | awk -v tag=GN -f scripts/awk/bam_clearTag.awk - \
            | awk -v tag=GX -f scripts/awk/bam_clearTag.awk - \
            | awk -v N={params.TRIM_N} -v option={params.TRIM_OPTION} -f scripts/awk/bam_trimNBases.awk - \
            | {EXEC['SAMTOOLS']} view -bS \
            > {output.BAM}
            """
        )


# Run bowtie2 on miRNA reference from pirbase in end-to-end/sensitive mode
#   https://bowtie-bio.sourceforge.net/bowtie2/manual.shtml#end-to-end-alignment-versus-local-alignment
#   change to `--sensitive-local` for other applications
# To generate: `bowtie2-build mmu.gold.fa.gz ./index > build.log`
rule bowtie2_align_mature_miRNA:
    input:
        BAM="{OUTDIR}/{SAMPLE}/miRNA/{RECIPE}/tmp.bam",
        # R1_FQ_FILTERED = '{OUTDIR}/{SAMPLE}/tmp/{SAMPLE}_R1_final_filtered_short.fq.gz',
        # R2_FQ_FILTERED = '{OUTDIR}/{SAMPLE}/tmp/{SAMPLE}_R2_final_filtered_short.fq.gz'        
    output:
        TMP_ALIGNED_BAM=temp("{OUTDIR}/{SAMPLE}/miRNA/{RECIPE}/tmp.aligned.bam"),
        MATURE_BAM="{OUTDIR}/{SAMPLE}/miRNA/{RECIPE}/mature.aligned.bam",
        UNALIGNED_BAM=temp("{OUTDIR}/{SAMPLE}/miRNA/{RECIPE}/unaligned.bam"),
    params:
        OUTDIR=config["OUTDIR"],
        MEMLIMIT=config["MEMLIMIT"],
        REF=config["miRNA_MATURE_INDEX"],
    log:
        log="{OUTDIR}/{SAMPLE}/miRNA/{RECIPE}/bowtie2_mature.log",
    # 1
    threads: config["CORES"]
    run:
        # Align to mature miR reference
        shell(
            f"""
            {EXEC['BOWTIE2']} \
                -x {params.REF} \
                -b {input.BAM} \
                -p {threads} \
                --preserve-tags \
            2> {log.log} \
            | {EXEC['SAMTOOLS']} view -bS \
            > {output.TMP_ALIGNED_BAM}
            """
        )
        # --very-sensitive-local \

        # Filter to save aligned reads (mature miRs)
        shell(
            f"""
            {EXEC['SAMTOOLS']} view \
                -b -F 4 \
                {output.TMP_ALIGNED_BAM} \
                > {output.MATURE_BAM}
            """
        )

        # Get unaligned reads for hairpin alignment
        shell(
            f"""
            {EXEC['SAMTOOLS']} view -f 4 {output.TMP_ALIGNED_BAM} \
            | awk -f scripts/awk/bam_clearAlignment.awk - \
            | awk -v tag=AS -f scripts/awk/bam_clearTag.awk - \
            | awk -v tag=XS -f scripts/awk/bam_clearTag.awk - \
            | {EXEC['SAMTOOLS']} view -bS \
            > {output.UNALIGNED_BAM}
            """
        )


# Align to hairpin miR reference, toss unaligned
rule bowtie2_align_hairpin_miRNA:
    input:
        BAM="{OUTDIR}/{SAMPLE}/miRNA/{RECIPE}/unaligned.bam",
    output:
        BAM="{OUTDIR}/{SAMPLE}/miRNA/{RECIPE}/hairpin.aligned.bam",
    params:
        OUTDIR=config["OUTDIR"],
        MEMLIMIT=config["MEMLIMIT"],
        REF=config["miRNA_HAIRPIN_INDEX"],
    log:
        log="{OUTDIR}/{SAMPLE}/miRNA/{RECIPE}/bowtie2_hairpin.log",
    # 1
    threads: config["CORES"]
    run:
        shell(
            f"""
            {EXEC['BOWTIE2']} \
                -x {params.REF} \
                -b {input.BAM} \
                -p {threads} \
                --no-unal \
                --preserve-tags \
            2> {log.log} \
            | {EXEC['SAMTOOLS']} view -bS \
            > {output.BAM}
            """
        )
        # --very-sensitive-local \



# Merge hairpin & mature alignment records
rule merge_aligned_bams_miRNA:
    input:
        MATURE_BAM="{OUTDIR}/{SAMPLE}/miRNA/{RECIPE}/mature.aligned.bam",
        HAIRPIN_BAM="{OUTDIR}/{SAMPLE}/miRNA/{RECIPE}/hairpin.aligned.bam",
    output:
        BAM="{OUTDIR}/{SAMPLE}/miRNA/{RECIPE}/aligned.bam",
    params:
        OUTDIR=config["OUTDIR"],
        MEMLIMIT=config["MEMLIMIT"],
        REF=config["miRNA_HAIRPIN_INDEX"],
    threads: 1
    # config['CORES']
    run:
        shell(
            f"""
            {EXEC['SAMTOOLS']} merge {output.BAM} {input.HAIRPIN_BAM} {input.MATURE_BAM}
            """
        )


# Index the deduplicated .bam file
rule sortAlignedBAM_miRNA:
    input:
        BAM="{OUTDIR}/{SAMPLE}/miRNA/{RECIPE}/aligned.bam",
    output:
        BAM=temp("{OUTDIR}/{SAMPLE}/miRNA/{RECIPE}/aligned.sorted.bam"),
    threads: config["CORES"]
    run:
        shell(
            f"""
            {EXEC['SAMTOOLS']} sort -@ {threads} {input.BAM} \
            > {output.BAM}
            """
        )


# Tag bam w/ chromosome/miRNA it aligned to (saved in "GN" tag, like STARsolo)
rule tagSortedBam_miRNA:
    input:
        BAM="{OUTDIR}/{SAMPLE}/miRNA/{RECIPE}/aligned.sorted.bam",
    output:
        BAM="{OUTDIR}/{SAMPLE}/miRNA/{RECIPE}/aligned.sorted.tagged.bam",  #TODO: add temp() in favor of just keeping the deduped bam?
    params:
        OUTDIR=config["OUTDIR"],
    threads: 2
    run:
        shell(
            f"""
            {EXEC['SAMTOOLS']} view -h {input.BAM} \
            | awk -f scripts/awk/bam_chr2tag.awk -v tag=GN - \
            | {EXEC['SAMTOOLS']} view -bS - \
            > {output.BAM}
            """
        )


# Index the sorted & deduplicated .bam file
rule indexSortedTaggedBAM_miRNA:
    input:
        BAM="{OUTDIR}/{SAMPLE}/miRNA/{RECIPE}/aligned.sorted.tagged.bam",
    output:
        BAI=temp("{OUTDIR}/{SAMPLE}/miRNA/{RECIPE}/aligned.sorted.tagged.bam.bai"),
    threads: config["CORES"]
    run:
        shell(
            f"""
            {EXEC['SAMTOOLS']} index -@ {threads} {input.BAM}
            """
        )


# Generate count matrix w/ umi-tools for miRNAs
rule umitools_count_miRNA:
    input:
        BAM="{OUTDIR}/{SAMPLE}/miRNA/{RECIPE}/aligned.sorted.tagged.bam",
        BAI="{OUTDIR}/{SAMPLE}/miRNA/{RECIPE}/aligned.sorted.tagged.bam.bai",
    output:
        COUNTS="{OUTDIR}/{SAMPLE}/miRNA/{RECIPE}/counts.tsv.gz",
    params:
        OUTDIR=config["OUTDIR"],
    log:
        log="{OUTDIR}/{SAMPLE}/miRNA/{RECIPE}/count.log",
    threads: 1
    run:
        shell(
            f"""
            {EXEC['UMITOOLS']} count \
                --extract-umi-method=tag \
                --per-gene \
                --per-cell \
                --cell-tag=CB \
                --gene-tag=GN \
                --umi-tag=UB \
                --log={log.log} \
                -I {input.BAM} \
                -S {output.COUNTS}
            """
        )
        # --wide-format-cell-counts \



# Convert the long-format counts into a format that people can actually use
rule counts_to_sparse_miRNA:
    input:
        COUNTS="{OUTDIR}/{SAMPLE}/miRNA/{RECIPE}/counts.tsv.gz",
    output:
        BCS="{OUTDIR}/{SAMPLE}/miRNA/{RECIPE}/raw/barcodes.tsv.gz",
        FEATS="{OUTDIR}/{SAMPLE}/miRNA/{RECIPE}/raw/features.tsv.gz",
        COUNTS="{OUTDIR}/{SAMPLE}/miRNA/{RECIPE}/raw/matrix.mtx.gz",
    params:
        OUTDIR=config["OUTDIR"],
    threads: 1
    run:
        output_dir = output.COUNTS.replace("/matrix.mtx.gz", "")
        shell(
            f"""
            mkdir -p {output_dir}
            python scripts/py/long2mtx.py {input.COUNTS} {output_dir}
            """
        )


# Dedup the .bam (do NOT split across chromosomes, b/c of custom reference)
rule umitools_dedupSortedBAM_miRNA:
    input:
        BB_WHITELIST="{OUTDIR}/{SAMPLE}/bb/whitelist.txt",
        BAM="{OUTDIR}/{SAMPLE}/miRNA/{RECIPE}/aligned.sorted.tagged.bam",
    output:
        BAM="{OUTDIR}/{SAMPLE}/miRNA/{RECIPE}/aligned.sorted.tagged.dedup.bam",
    threads: config["CORES"]
    log:
        log="{OUTDIR}/{SAMPLE}/miRNA/{RECIPE}/dedup.log",
    run:
        tmp_recipe = RECIPE_DICT[wildcards.SAMPLE]

        whitelist = input.BB_WHITELIST

        shell(
            f"""
            bash scripts/bash/dedup.sh \
                {input.BAM} \
                {whitelist} \
                {threads} \
                {output.BAM} \
                {OUTDIR}/{wildcards.SAMPLE}/tmp/dedup \
                {log.log}
            """
        )


# Index the sorted & deduplicated .bam file
rule indexSortedTaggedDedupBAM_miRNA:
    input:
        BAM="{OUTDIR}/{SAMPLE}/miRNA/{RECIPE}/aligned.sorted.tagged.dedup.bam",
    output:
        BAI="{OUTDIR}/{SAMPLE}/miRNA/{RECIPE}/aligned.sorted.tagged.dedup.bam.bai",
    threads: config["CORES"]
    run:
        shell(
            f"""
            {EXEC['SAMTOOLS']} index -@ {threads} {input.BAM}
            """
        )
