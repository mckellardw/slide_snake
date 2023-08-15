# Remove reads that don't have a corrected spot/cell barcode with samtools, then remove duplicates w/ **umi-tools**
## High mem usage? Check here! https://umi-tools.readthedocs.io/en/latest/faq.html
#TODO: dedup after strand-splitting? 
#TODO: add exec paths for samtools, bamtools
rule umitools_dedupBAM:
    input:
        BB_WHITELIST = "{OUTDIR}/{sample}/bb/whitelist.txt",
        BB_1 = "{OUTDIR}/{sample}/bb/whitelist_1.txt",
        BB_2 = "{OUTDIR}/{sample}/bb/whitelist_2.txt",
        SORTEDBAM = '{OUTDIR}/{sample}/STARsolo/Aligned.sortedByCoord.out.bam'
    output:
        DEDUPBAM = '{OUTDIR}/{sample}/STARsolo/Aligned.sortedByCoord.dedup.out.bam'
    threads:
        config['CORES']
    log:
        '{OUTDIR}/{sample}/dedup.log'
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

# Index the deduplicated .bam file
rule umitools_indexDedupBAM:
    input:
        SORTEDBAM = '{OUTDIR}/{sample}/Aligned.sortedByCoord.dedup.out.bam'
    output:
        BAI = '{OUTDIR}/{sample}/Aligned.sortedByCoord.dedup.out.bam.bai'
    threads:
        config['CORES']
    run:
        shell(
            f"""
            {SAMTOOLS_EXEC} index -@ {threads} {input.SORTEDBAM}
            """
        )

# Split .bam file by strand for IGV browsing
rule strand_split_dedup_bam:
    input:
        DEDUPBAM = '{OUTDIR}/{sample}/Aligned.sortedByCoord.dedup.out.bam'
    output:
        FWDBAM = '{OUTDIR}/{sample}/Aligned.sortedByCoord.dedup.out.fwd.bam',
        REVBAM = '{OUTDIR}/{sample}/Aligned.sortedByCoord.dedup.out.rev.bam'
    threads:
        1
    run:
        shell(
            f"""
            {SAMTOOLS_EXEC} view -b -F 0x10 {input.DEDUPBAM} > {output.FWDBAM}
            {SAMTOOLS_EXEC} view -b -f 0x10 {input.DEDUPBAM} > {output.REVBAM}
            """
        )

# Index the split/deduped bam files
rule indexSplitBAMs:
    input:
        FWDBAM = '{OUTDIR}/{sample}/Aligned.sortedByCoord.dedup.out.fwd.bam',
        REVBAM = '{OUTDIR}/{sample}/Aligned.sortedByCoord.dedup.out.rev.bam'
    output:
        FWDBAI = '{OUTDIR}/{sample}/Aligned.sortedByCoord.dedup.out.fwd.bam.bai',
        REVBAI = '{OUTDIR}/{sample}/Aligned.sortedByCoord.dedup.out.rev.bam.bai'
    threads:
        config['CORES']
    run:
        shell(
            f"""
            {SAMTOOLS_EXEC} index -@ {threads} {input.FWDBAM}
            {SAMTOOLS_EXEC} index -@ {threads} {input.REVBAM}
            """
        )