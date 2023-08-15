# bowtie2 documentation: https://bowtie-bio.sourceforge.net/bowtie2/manual.shtml

# pirbase download link - http://bigdata.ibp.ac.cn/piRBase/download.php
# *Note* - using gold standard piRNA refs ()

# Use STAR-barcoded .bam file & filter to remove reads longer than any piRNAs (##bp)
rule bowtie2_prep_bam:
    input:
        BAM = '{OUTDIR}/{sample}/STARsolo/Aligned.sortedByCoord.out.bam'
    output:
        BAM = temp('{OUTDIR}/{sample}/pirna/tmp.bam')
    params:
        OUTDIR = config['OUTDIR'],
        MEMLIMIT = config['MEMLIMIT'],
        MAX_SHORT_READ_LENGTH = config['MAX_SHORT_READ_LENGTH']
    threads:
        config['CORES']
    run:
        shell(
            f"""
            mkdir -p {OUTDIR}/{wildcards.sample}/pirna

            {SAMTOOLS_EXEC} view {input.bam} \
            | awk 'length($10) > {params.MAX_SHORT_READ_LENGTH} || $1 ~ /^@/' \
            | awk 'BEGIN {{FS=OFS=\"\\t\"}} !/^@/ {{${3}=\"*\"; ${4}=\"0\"; ${5}=\"0\";t \$6=\"*\"; ${7}=\"*\"; ${8}=\"0\"; ${9}=\"0\"}} {{print}}' \
            | {SAMTOOLS_EXEC} view -bS {output.BAM}
            """
        )
# samtools view -h input.bam | awk 'length(\$10) > 30 || \$1 ~ /^@/' | samtools view -bS - > output.bam

# samtools view -h aligned.bam | awk 'BEGIN {FS=OFS="\t"} !/^@/ {\$3="*"; \$4="0"; \$5="0"; \$6="*"; \$7="*"; \$8="0"; \$9="0"} {print}' | samtools view -b -o unaligned.bam -

# Run bowtie2 on piRNA reference from pirbase
# To generate: `bowtie2-build mmu.gold.fa.gz ./index > build.log`
rule bowtie2_piRNA:
    input:
        BAM = '{OUTDIR}/{sample}/pirna/tmp.bam'
        # R1_FQ_FILTERED = '{OUTDIR}/{sample}/tmp/{sample}_R1_final_filtered_short.fq.gz',
        # R2_FQ_FILTERED = '{OUTDIR}/{sample}/tmp/{sample}_R2_final_filtered_short.fq.gz'        
    output:
        BAM = '{OUTDIR}/{sample}/pirna/aligned.bam'
    params:
        OUTDIR = config['OUTDIR'],
        MEMLIMIT = config['MEMLIMIT'],
        REF = '/workdir/dwm269/genomes/pirbase/mmu/gold/bowtie2' # just mouse implemented now...
    log:
        '{OUTDIR}/{sample}/pirna/bowtie2.log'    
    threads:
        config['CORES']
    run:
        shell(
            f"""
            {BOWTIE2_EXEC} \
            -x {params.REF} \
            -b {input.bam_file} \
            -p {threads} \
            --sensitive-local \
            --preserve-tags \
            --no-unal \
            --met-file {log} \
            | {SAMTOOLS_EXEC} view -bS - > {output}
            """
        )


# Index the deduplicated .bam file
rule umitools_sortDedupBAM_piRNA:
    input:
        BAM = '{OUTDIR}/{sample}/pirna/aligned.bam'
    output:
        BAM = '{OUTDIR}/{sample}/pirna/aligned.sorted.bam'
    threads:
        config['CORES']
    run:
        shell(
            f"""
            {SAMTOOLS_EXEC} sort -@ {threads} {input.BAM} > {output.BAM}
            """
        )


# Dedup the bam
rule umitools_dedupSortedBAM_piRNA:
    input:
        BB_WHITELIST = "{OUTDIR}/{sample}/bb/whitelist.txt",
        BAM = '{OUTDIR}/{sample}/pirna/aligned.sorted.bam'
    output:
        BAM = '{OUTDIR}/{sample}/pirna/aligned.sorted.dedup.bam'
    threads:
        config['CORES']
    log:
        '{OUTDIR}/{sample}/pirna/dedup.log'
    run:
        tmp_recipe = RECIPE_DICT[wildcards.sample]

        whitelist = input.BB_WHITELIST

        shell(
            f"""
            bash scripts/split_dedup.sh \
            {input.SORTEDBAM} \
            {whitelist} \
            {threads} \
            {output.DEDUPBAM} \
            {OUTDIR}/{wildcards.sample}/tmp/dedup \
            | tee {log}
            """
        )


# Index the sorted & deduplicated .bam file
rule umitools_indexSortedDedupBAM_piRNA:
    input:
        BAM = '{OUTDIR}/{sample}/pirna/aligned.sorted.dedup.bam'
    output:
        BAI = '{OUTDIR}/{sample}/pirna/aligned.sorted.dedup.bam.bai'
    threads:
        config['CORES']
    run:
        shell(
            f"""
            {SAMTOOLS_EXEC} index -@ {threads} {input.BAM}
            """
        )

# Tag bam w/ chromosome/piRNA it aligned to, 
#   then generate count matrix w/ umi-tools for piRNAs
rule bowtie2_piRNA_counts:
    input:
        BAM = '{OUTDIR}/{sample}/pirna/aligned.sorted.dedup.bam'
    output:
        COUNTS = '{OUTDIR}/{sample}/pirna/counts.tsv.gz'
    params:
        OUTDIR = config['OUTDIR']
    log:
        '{OUTDIR}/{sample}/pirna/count.log'
    run:
        shell(
            f"""
            {SAMTOOLS_EXEC} view {input.BAM} \
            | awk -f scripts/add_bt_tag.awk \
            | {SAMTOOLS_EXEC} view -b \
            | {UMITOOLS_EXEC} count \
            --per-gene \
            --per-cell \
            --cell-barcode-tag=CB \
            --gene-tag=BT \
            --log={log} \
            -I {input.BAM} \
            -S {output.COUNTS}
            """
        )