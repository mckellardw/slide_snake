####################################################
# Summary of Workflow:

# 1) bowtie2_prep_bam: Filter and preprocess the STAR-aligned BAM file to remove 
#       reads longer than piRNAs, and clear alignment information while keeping tags.

# 2) bowtie2_piRNA: Align the preprocessed BAM file to the piRNA reference using 
#       bowtie2, and save the aligned BAM file.

# 3) umitools_sortAlignedBAM_piRNA: Sort the aligned BAM file using samtools.

# 4) umitools_dedupSortedBAM_piRNA: Deduplicate the sorted BAM file using umi-tools.

# 5) umitools_indexSortedDedupBAM_piRNA: Index the deduplicated BAM file using samtools.

# 6) bowtie2_piRNA_counts: Tag the deduplicated BAM file with chromosome/piRNA information 
#       and generate a count matrix using umi-tools.


# bowtie2 documentation: https://bowtie-bio.sourceforge.net/bowtie2/manual.shtml
# pirbase download link - http://bigdata.ibp.ac.cn/piRBase/download.php
#       *Note* - using gold standard piRNA refs ()
####################################################

# Use STAR-barcoded .bam file & filter:
#   - Remove reads longer than any piRNAs (34bp)
#   - Remove reads missing a cell barcode ("CB" tag) or UMI ("UB" tag)
#   - Remove all aligment info, but keep tags
rule bowtie2_prep_bam_piRNA:
    input:
        BAM = '{OUTDIR}/{sample}/STARsolo/Aligned.sortedByCoord.out.bam'
    output:
        BAM = temp('{OUTDIR}/{sample}/piRNA/tmp.bam')
    params:
        OUTDIR = config['OUTDIR'],
        MEMLIMIT = config['MEMLIMIT'],
        MAX_SHORT_READ_LENGTH = config['MAX_SHORT_READ_LENGTH']
    threads:
        1
    run:
        shell(
            f"""
            mkdir -p {OUTDIR}/{wildcards.sample}/piRNA

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

# Run bowtie2 on piRNA reference from pirbase in end-to-end/sensitive mode 
#   https://bowtie-bio.sourceforge.net/bowtie2/manual.shtml#end-to-end-alignment-versus-local-alignment
#   change to `--sensitive-local` for other applications
# To generate: `bowtie2-build mmu.gold.fa.gz ./index > build.log`
rule bowtie2_align_piRNA:
    input:
        BAM = '{OUTDIR}/{sample}/piRNA/tmp.bam'
        # R1_FQ_FILTERED = '{OUTDIR}/{sample}/tmp/{sample}_R1_final_filtered_short.fq.gz',
        # R2_FQ_FILTERED = '{OUTDIR}/{sample}/tmp/{sample}_R2_final_filtered_short.fq.gz'        
    output:
        BAM = '{OUTDIR}/{sample}/piRNA/aligned.bam'
    params:
        OUTDIR = config['OUTDIR'],
        MEMLIMIT = config['MEMLIMIT'],
        REF = config['piRNA_INDEX']
    log:
        '{OUTDIR}/{sample}/piRNA/bowtie2.log'    
    threads:
        # 1
        config['CORES']
    run:
        shell(
            f"""
            {BOWTIE2_EXEC} \
            -x {params.REF} \
            -b {input.BAM} \
            -p {threads} \
            --very-sensitive-local \
            --preserve-tags \
            --no-unal \
            2> {log} \
            | {SAMTOOLS_EXEC} view -bS > {output.BAM}
            """
        )


# Index the deduplicated .bam file
rule sortAlignedBAM_piRNA:
    input:
        BAM = '{OUTDIR}/{sample}/piRNA/aligned.bam'
    output:
        BAM = temp('{OUTDIR}/{sample}/piRNA/aligned.sorted.bam')
    threads:
        config['CORES']
    run:
        shell(
            f"""
            {SAMTOOLS_EXEC} sort -@ {threads} {input.BAM} > {output.BAM}
            """
        )


# Tag bam w/ chromosome/piRNA it aligned to
rule tagSortedBam_piRNA:
    input:
        BAM = '{OUTDIR}/{sample}/piRNA/aligned.sorted.bam'
    output:
        BAM = '{OUTDIR}/{sample}/piRNA/aligned.sorted.tagged.bam' #TODO: add temp() in favor of just keeping the deduped bam?
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
rule indexSortedTaggedBAM_piRNA:
    input:
        BAM = '{OUTDIR}/{sample}/piRNA/aligned.sorted.tagged.bam'
    output:
        BAI = temp('{OUTDIR}/{sample}/piRNA/aligned.sorted.tagged.bam.bai')
    threads:
        config['CORES']
    run:
        shell(
            f"""
            {SAMTOOLS_EXEC} index -@ {threads} {input.BAM}
            """
        )

# Generate count matrix w/ umi-tools for piRNAs
rule umitools_count_piRNA:
    input:
        BAM = '{OUTDIR}/{sample}/piRNA/aligned.sorted.tagged.bam',
        BAI = '{OUTDIR}/{sample}/piRNA/aligned.sorted.tagged.bam.bai'
    output:        
        COUNTS = '{OUTDIR}/{sample}/piRNA/counts.tsv.gz'
    params:
        OUTDIR = config['OUTDIR']
    log:
        '{OUTDIR}/{sample}/piRNA/count.log'
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
rule umitools_dedupSortedBAM_piRNA:
    input:
        BB_WHITELIST = "{OUTDIR}/{sample}/bb/whitelist.txt",
        BAM = '{OUTDIR}/{sample}/piRNA/aligned.sorted.tagged.bam'
    output:
        BAM = '{OUTDIR}/{sample}/piRNA/aligned.sorted.tagged.dedup.bam'
    threads:
        config['CORES']
    log:
        '{OUTDIR}/{sample}/piRNA/dedup.log'
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
rule indexSortedTaggedDedupBAM_piRNA:
    input:
        BAM = '{OUTDIR}/{sample}/piRNA/aligned.sorted.tagged.dedup.bam'
    output:
        BAI = '{OUTDIR}/{sample}/piRNA/aligned.sorted.tagged.dedup.bam.bai'
    threads:
        config['CORES']
    run:
        shell(
            f"""
            {SAMTOOLS_EXEC} index -@ {threads} {input.BAM}
            """
        )
