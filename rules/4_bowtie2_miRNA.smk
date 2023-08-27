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
#   - Remove reads longer than any miRNAs (34bp)
#   - Remove reads missing a cell barcode ("CB" tag) or UMI ("UB" tag)
#   - Remove all aligment info, but keep tags
rule bowtie2_prep_bam_miRNA:
    input:
        BAM = '{OUTDIR}/{sample}/STARsolo/Aligned.sortedByCoord.out.bam'
    output:
        BAM = temp('{OUTDIR}/{sample}/miRNA/tmp.bam')
    params:
        OUTDIR = config['OUTDIR'],
        MEMLIMIT = config['MEMLIMIT'],
        MAX_SHORT_READ_LENGTH = config['MAX_SHORT_READ_LENGTH']
    threads:
        1
    run:
        shell(
            f"""
            mkdir -p {OUTDIR}/{wildcards.sample}/miRNA

            {SAMTOOLS_EXEC} view {input.BAM} \
            | awk -f scripts/bam_shortPassReadFilter.awk -v max_length={params.MAX_SHORT_READ_LENGTH} \
            | awk -v tag=CB -f scripts/bam_filterEmptyTag.awk - \
            | awk -v tag=UB -f scripts/bam_filterEmptyTag.awk - \
            | awk -f scripts/bam_clearAlignment.awk - \
            | {SAMTOOLS_EXEC} view -bS > {output.BAM}
            """
        )
# samtools view -h input.bam | awk 'length(\$10) > 30 || \$1 ~ /^@/' | samtools view -bS - > output.bam
# samtools view -h aligned.bam | awk 'BEGIN {FS=OFS="\t"} !/^@/ {\$3="*"; \$4="0"; \$5="0"; \$6="*"; \$7="*"; \$8="0"; \$9="0"} {print}' | samtools view -b -o unaligned.bam -

# Run bowtie2 on miRNA reference from pirbase in end-to-end/sensitive mode 
#   https://bowtie-bio.sourceforge.net/bowtie2/manual.shtml#end-to-end-alignment-versus-local-alignment
#   change to `--sensitive-local` for other applications
# To generate: `bowtie2-build mmu.gold.fa.gz ./index > build.log`
rule bowtie2_align_mature_miRNA:
    input:
        BAM = '{OUTDIR}/{sample}/miRNA/tmp.bam'
        # R1_FQ_FILTERED = '{OUTDIR}/{sample}/tmp/{sample}_R1_final_filtered_short.fq.gz',
        # R2_FQ_FILTERED = '{OUTDIR}/{sample}/tmp/{sample}_R2_final_filtered_short.fq.gz'        
    output:
        TMP_ALIGNED_BAM = temp('{OUTDIR}/{sample}/miRNA/tmp.aligned.bam'),
        MATURE_BAM = '{OUTDIR}/{sample}/miRNA/mature.aligned.bam',
        UNALIGNED_BAM = temp('{OUTDIR}/{sample}/miRNA/unaligned.bam')
    params:
        OUTDIR = config['OUTDIR'],
        MEMLIMIT = config['MEMLIMIT'],
        REF = config['miRNA_MATURE_INDEX']
    log:
        '{OUTDIR}/{sample}/miRNA/bowtie2_mature.log'    
    threads:
        1
        # config['CORES']
    run:
        # Align to mature miR reference
        shell(
            f"""
            {BOWTIE2_EXEC} \
            -x {params.REF} \
            -b {input.BAM} \
            -p {threads} \
            --very-sensitive-local \
            --preserve-tags \
            2> {log} \
            | {SAMTOOLS_EXEC} view -bS \
            > {output.TMP_ALIGNED_BAM}
            """
        )

        # Filter to save aligned reads (mature miRs)
        shell(
            f"""
            {SAMTOOLS_EXEC} view -b -F 4 {output.TMP_ALIGNED_BAM} > {output.MATURE_BAM}
            """
        )

        # Get unaligned reads for hairpin alignment
        shell(
            f"""
            {SAMTOOLS_EXEC} view -f 4 {output.TMP_ALIGNED_BAM} \
            | awk -f scripts/bam_clearAlignment.awk - \
            | {SAMTOOLS_EXEC} view -bS \
            > {output.UNALIGNED_BAM}
            """
        )

# Align to hairpin miR reference, toss unaligned
rule bowtie2_align_hairpin_miRNA:
    input:
        BAM = '{OUTDIR}/{sample}/miRNA/unaligned.bam'      
    output:
        BAM = '{OUTDIR}/{sample}/miRNA/hairpin.aligned.bam'
    params:
        OUTDIR = config['OUTDIR'],
        MEMLIMIT = config['MEMLIMIT'],
        REF = config['miRNA_HAIRPIN_INDEX']
    log:
        '{OUTDIR}/{sample}/miRNA/bowtie2_hairpin.log'    
    threads:
        1
        # config['CORES']
    run:
        shell(
            f"""
            {BOWTIE2_EXEC} \
            -x {params.REF} \
            -b {input.BAM} \
            -p {threads} \
            --very-sensitive-local \
            --no-unal \
            --preserve-tags \
            2> {log} \
            | {SAMTOOLS_EXEC} view -bS > {output.BAM}
            """
        )

# Merge hairpin & mature alignment records
rule merge_aligned_bams_miRNA:
    input:
        MATURE_BAM = '{OUTDIR}/{sample}/miRNA/mature.aligned.bam',
        HAIRPIN_BAM = '{OUTDIR}/{sample}/miRNA/hairpin.aligned.bam'
    output:
        BAM = '{OUTDIR}/{sample}/miRNA/aligned.bam'
    params:
        OUTDIR = config['OUTDIR'],
        MEMLIMIT = config['MEMLIMIT'],
        REF = config['miRNA_HAIRPIN_INDEX']
    log:
        '{OUTDIR}/{sample}/miRNA/bowtie2_hairpin.log'    
    threads:
        1
        # config['CORES']
    run:
        shell(
            f"""
            {SAMTOOLS_EXEC} merge {output.BAM} {input.HAIRPIN_BAM} {input.MATURE_BAM}
            """
        )


# Index the deduplicated .bam file
rule sortAlignedBAM_miRNA:
    input:
        BAM = '{OUTDIR}/{sample}/miRNA/aligned.bam'
    output:
        BAM = temp('{OUTDIR}/{sample}/miRNA/aligned.sorted.bam')
    threads:
        config['CORES']
    run:
        shell(
            f"""
            {SAMTOOLS_EXEC} sort -@ {threads} {input.BAM} > {output.BAM}
            """
        )


# Tag bam w/ chromosome/miRNA it aligned to
rule tagSortedBam_miRNA:
    input:
        BAM = '{OUTDIR}/{sample}/miRNA/aligned.sorted.bam'
    output:
        BAM = '{OUTDIR}/{sample}/miRNA/aligned.sorted.tagged.bam' #TODO: add temp() in favor of just keeping the deduped bam?
    params:
        OUTDIR = config['OUTDIR']
    threads:
        1
    run:
        shell(
            f"""
            {SAMTOOLS_EXEC} view -h {input.BAM} \
            | awk -f scripts/bam_chr2tag.awk - \
            | {SAMTOOLS_EXEC} view -bS - \
            > {output.BAM}
            """
        )


# Index the sorted & deduplicated .bam file
rule indexSortedTaggedBAM_miRNA:
    input:
        BAM = '{OUTDIR}/{sample}/miRNA/aligned.sorted.tagged.bam'
    output:
        BAI = temp('{OUTDIR}/{sample}/miRNA/aligned.sorted.tagged.bam.bai')
    threads:
        config['CORES']
    run:
        shell(
            f"""
            {SAMTOOLS_EXEC} index -@ {threads} {input.BAM}
            """
        )

# Generate count matrix w/ umi-tools for miRNAs
rule umitools_count_miRNA:
    input:
        BAM = '{OUTDIR}/{sample}/miRNA/aligned.sorted.tagged.bam',
        BAI = '{OUTDIR}/{sample}/miRNA/aligned.sorted.tagged.bam.bai'
    output:        
        COUNTS = '{OUTDIR}/{sample}/miRNA/counts.tsv.gz'
    params:
        OUTDIR = config['OUTDIR']
    log:
        '{OUTDIR}/{sample}/miRNA/count.log'
    run:
        shell(
            f"""
            {UMITOOLS_EXEC} count \
            --extract-umi-method=tag \
            --per-gene \
            --per-cell \
            --wide-format-cell-counts \
            --cell-tag=CB \
            --gene-tag=BT \
            --umi-tag=UB \
            --log={log} \
            -I {input.BAM} \
            -S {output.COUNTS}
            """
        )


# Dedup the .bam (do NOT split across chromosomes, b/c of custom reference)
rule umitools_dedupSortedBAM_miRNA:
    input:
        BB_WHITELIST = "{OUTDIR}/{sample}/bb/whitelist.txt",
        BAM = '{OUTDIR}/{sample}/miRNA/aligned.sorted.tagged.bam'
    output:
        BAM = '{OUTDIR}/{sample}/miRNA/aligned.sorted.tagged.dedup.bam'
    threads:
        config['CORES']
    log:
        '{OUTDIR}/{sample}/miRNA/dedup.log'
    run:
        tmp_recipe = RECIPE_DICT[wildcards.sample]

        whitelist = input.BB_WHITELIST

        shell(
            f"""
            bash scripts/dedup.sh \
            {input.BAM} \
            {whitelist} \
            {threads} \
            {output.BAM} \
            {OUTDIR}/{wildcards.sample}/tmp/dedup \
            {log}
            """
        )


# Index the sorted & deduplicated .bam file
rule indexSortedTaggedDedupBAM_miRNA:
    input:
        BAM = '{OUTDIR}/{sample}/miRNA/aligned.sorted.tagged.dedup.bam'
    output:
        BAI = '{OUTDIR}/{sample}/miRNA/aligned.sorted.tagged.dedup.bam.bai'
    threads:
        config['CORES']
    run:
        shell(
            f"""
            {SAMTOOLS_EXEC} index -@ {threads} {input.BAM}
            """
        )
