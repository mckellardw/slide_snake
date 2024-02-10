rule BLAST_align_mature_miRNA:
    input:
        BAM = '{OUTDIR}/{SAMPLE}/miRNA/tmp.bam'
        # R1_FQ_FILTERED = '{OUTDIR}/{SAMPLE}/tmp/{SAMPLE}_R1_final_filtered_short.fq.gz',
        # R2_FQ_FILTERED = '{OUTDIR}/{SAMPLE}/tmp/{SAMPLE}_R2_final_filtered_short.fq.gz'        
    output:
        TMP_ALIGNED_BAM = temp('{OUTDIR}/{SAMPLE}/miRNA/tmp.aligned.bam'),
        MATURE_BAM = '{OUTDIR}/{SAMPLE}/miRNA/mature.aligned.bam',
        UNALIGNED_BAM = temp('{OUTDIR}/{SAMPLE}/miRNA/unaligned.bam')
    params:
        OUTDIR = config['OUTDIR'],
        MEMLIMIT = config['MEMLIMIT'],
        REF = config['miRNA_MATURE_INDEX']
    log:
        log = '{OUTDIR}/{SAMPLE}/miRNA/BLAST_mature.log'    
    threads:
        1
        # config['CORES']
    run:

        # Convert BAM to FASTA
        shell(
            f"""
            {EXEC['SAMTOOLS']} fasta input.bam > input.fasta
            """
        )

        # Run BLAST
        shell(
            f"""
            {EXEC['BLASTN']} \
                -query input.fasta \
                -db your_db \
                -outfmt 6 \
                -out blast_results.txt
            """
        )


        # Parse BLAST results and add to BAM file
        shell(
            f"""
            awk 'BEGIN {OFS="\t"} {print $1, "BL:Z:" $2}' blast_results.txt > blast_tags.txt
            """
        )

        # Add tags to BAM file
        shell(
            f"""
            {EXEC['SAMTOOLS']} view -h input.bam \
            | awk -v tags_file="blast_tags.txt" 'BEGIN {while ((getline < tags_file) > 0) tags[$1]=$2} {if ($1 in tags) print $0 "\t" tags[$1]; else print $0}' \
            | {EXEC['SAMTOOLS']} view -b \
            > output.bam
            """
        )

        # Filter to save aligned reads (mature miRs)
        shell(
            f"""
            {EXEC['SAMTOOLS']} view -b -F 4 {output.TMP_ALIGNED_BAM} > {output.MATURE_BAM}
            """
        )

        # Get unaligned reads for hairpin alignment
        shell(
            f"""
            {EXEC['SAMTOOLS']} view -f 4 {output.TMP_ALIGNED_BAM} \
            | awk -f scripts/awk/bam_clearAlignment.awk - \
            | {EXEC['SAMTOOLS']} view -bS \
            > {output.UNALIGNED_BAM}
            """
        )

# Align to hairpin miR reference, toss unaligned
rule BLAST_align_hairpin_miRNA:
    input:
        BAM = '{OUTDIR}/{SAMPLE}/miRNA/unaligned.bam'      
    output:
        BAM = '{OUTDIR}/{SAMPLE}/miRNA/hairpin.aligned.bam'
    params:
        OUTDIR = config['OUTDIR'],
        MEMLIMIT = config['MEMLIMIT'],
        REF = config['miRNA_HAIRPIN_INDEX']
    log:
        log = '{OUTDIR}/{SAMPLE}/miRNA/BLAST_hairpin.log'    
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
            2> {log.log} \
            | {EXEC['SAMTOOLS']} view -bS > {output.BAM}
            """
        )

# Merge hairpin & mature alignment records
rule merge_aligned_bams_miRNA:
    input:
        MATURE_BAM = '{OUTDIR}/{SAMPLE}/miRNA/mature.aligned.bam',
        HAIRPIN_BAM = '{OUTDIR}/{SAMPLE}/miRNA/hairpin.aligned.bam'
    output:
        BAM = '{OUTDIR}/{SAMPLE}/miRNA/aligned.bam'
    params:
        OUTDIR = config['OUTDIR'],
        MEMLIMIT = config['MEMLIMIT'],
        REF = config['miRNA_HAIRPIN_INDEX']
    log:
        '{OUTDIR}/{SAMPLE}/miRNA/BLAST_hairpin.log'    
    threads:
        1
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
        BAM = '{OUTDIR}/{SAMPLE}/miRNA/aligned.bam'
    output:
        BAM = temp('{OUTDIR}/{SAMPLE}/miRNA/aligned.sorted.bam')
    threads:
        config['CORES']
    run:
        shell(
            f"""
            {EXEC['SAMTOOLS']} sort -@ {threads} {input.BAM} > {output.BAM}
            """
        )


# Tag bam w/ chromosome/miRNA it aligned to
rule tagSortedBam_miRNA:
    input:
        BAM = '{OUTDIR}/{SAMPLE}/miRNA/aligned.sorted.bam'
    output:
        BAM = '{OUTDIR}/{SAMPLE}/miRNA/aligned.sorted.tagged.bam' #TODO: add temp() in favor of just keeping the deduped bam?
    params:
        OUTDIR = config['OUTDIR']
    threads:
        1
    run:
        shell(
            f"""
            {EXEC['SAMTOOLS']} view -h {input.BAM} \
            | awk -f scripts/awk/bam_chr2tag.awk - \
            | {EXEC['SAMTOOLS']} view -bS - \
            > {output.BAM}
            """
        )


# Index the sorted & deduplicated .bam file
rule indexSortedTaggedBAM_miRNA:
    input:
        BAM = '{OUTDIR}/{SAMPLE}/miRNA/aligned.sorted.tagged.bam'
    output:
        BAI = temp('{OUTDIR}/{SAMPLE}/miRNA/aligned.sorted.tagged.bam.bai')
    threads:
        config['CORES']
    run:
        shell(
            f"""
            {EXEC['SAMTOOLS']} index -@ {threads} {input.BAM}
            """
        )

# Generate count matrix w/ umi-tools for miRNAs
rule umitools_count_miRNA:
    input:
        BAM = '{OUTDIR}/{SAMPLE}/miRNA/aligned.sorted.tagged.bam',
        BAI = '{OUTDIR}/{SAMPLE}/miRNA/aligned.sorted.tagged.bam.bai'
    output:        
        COUNTS = '{OUTDIR}/{SAMPLE}/miRNA/counts.tsv.gz'
    params:
        OUTDIR = config['OUTDIR']
    log:
        log = '{OUTDIR}/{SAMPLE}/miRNA/count.log'
    run:
        shell(
            f"""
            {EXEC['UMITOOLS']} count \
                --extract-umi-method=tag \
                --per-gene \
                --per-cell \
                --wide-format-cell-counts \
                --cell-tag=CB \
                --gene-tag=BT \
                --umi-tag=UB \
                --log={log.log} \
                -I {input.BAM} \
                -S {output.COUNTS}
            """
        )


# Dedup the .bam (do NOT split across chromosomes, b/c of custom reference)
rule umitools_dedupSortedBAM_miRNA:
    input:
        BB_WHITELIST = "{OUTDIR}/{SAMPLE}/bb/whitelist.txt",
        BAM = '{OUTDIR}/{SAMPLE}/miRNA/aligned.sorted.tagged.bam'
    output:
        BAM = '{OUTDIR}/{SAMPLE}/miRNA/aligned.sorted.tagged.dedup.bam'
    threads:
        config['CORES']
    log:
        log = '{OUTDIR}/{SAMPLE}/miRNA/dedup.log'
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
        BAM = '{OUTDIR}/{SAMPLE}/miRNA/aligned.sorted.tagged.dedup.bam'
    output:
        BAI = '{OUTDIR}/{SAMPLE}/miRNA/aligned.sorted.tagged.dedup.bam.bai'
    threads:
        config['CORES']
    run:
        shell(
            f"""
            {EXEC['SAMTOOLS']} index -@ {threads} {input.BAM}
            """
        )
