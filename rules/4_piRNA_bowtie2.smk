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
        BAM = '{OUTDIR}/{sample}/STARsolo/Aligned.sortedByCoord.dedup.out.bam'
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
            | awk -v tag=AS -f scripts/bam_clearTag.awk - \
            | awk -v tag=GN -f scripts/bam_clearTag.awk - \
            | awk -v tag=GX -f scripts/bam_clearTag.awk - \
            | {SAMTOOLS_EXEC} view -bS > {output.BAM}
            """
        )


# Run bowtie2 on piRNA reference from pirbase in end-to-end/sensitive mode 
#   https://bowtie-bio.sourceforge.net/bowtie2/manual.shtml#end-to-end-alignment-versus-local-alignment
#   change to `--sensitive-local` for other applications
# To generate: `bowtie2-build mmu.gold.fa.gz ./index > build.log`
rule bowtie2_align_piRNA:
    input:
        BAM = '{OUTDIR}/{sample}/piRNA/tmp.bam'       
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
            | awk -f scripts/bam_chr2tag.awk -v tag=GN - \
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
            --cell-tag=CB \
            --gene-tag=GN \
            --umi-tag=UB \
            --log={log} \
            -I {input.BAM} \
            -S {output.COUNTS}
            """
        )

# Convert the long-format counts into a format that people can actually use
rule counts_to_sparse_piRNA:
    input:
        COUNTS = '{OUTDIR}/{sample}/piRNA/counts.tsv.gz'
    output:
        BCS = '{OUTDIR}/{sample}/piRNA/raw/barcodes.tsv.gz',
        FEATS = '{OUTDIR}/{sample}/piRNA/raw/features.tsv.gz',
        COUNTS = '{OUTDIR}/{sample}/piRNA/raw/matrix.mtx.gz'
    params:
        OUTDIR = config['OUTDIR']
    threads:
        1
    run:
        output_dir = output.COUNTS.replace('/matrix.mtx.gz','')
        shell(
            f"""
            mkdir -p {output_dir}
            python scripts/py/long2mtx.py {input.COUNTS} {output_dir}
            """   
        )


# Dedup the .bam (do NOT split across chromosomes, b/c of custom reference)
rule umitools_dedupSortedBAM_piRNA:
    input:
        BB_WHITELIST = '{OUTDIR}/{sample}/bb/whitelist.txt',
        BAM = '{OUTDIR}/{sample}/piRNA/aligned.sorted.tagged.bam'
    output:
        BAM = '{OUTDIR}/{sample}/piRNA/aligned.sorted.tagged.dedup.bam'
    threads:
        config['CORES']
    log:
        '{OUTDIR}/{sample}/piRNA/dedup.log'
    run:
        tmp_recipe = RECIPE_DICT[wildcards.sample]

        shell(
            f"""
            bash scripts/dedup.sh \
            {input.BAM} \
            {input.BB_WHITELIST} \
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
