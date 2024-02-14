# Remove reads that don't have a corrected spot/cell barcode with samtools, then remove duplicates w/ **umi-tools**
## High mem usage? Check here! https://umi-tools.readthedocs.io/en/latest/faq.html
#TODO: dedup after strand-splitting? 
#TODO: add exec paths for samtools, bamtools
rule umitools_dedupBAM:
    input:
        BB_WHITELIST = "{OUTDIR}/{SAMPLE}/bb/whitelist.txt",
        BB_1 = "{OUTDIR}/{SAMPLE}/bb/whitelist_1.txt",
        BB_2 = "{OUTDIR}/{SAMPLE}/bb/whitelist_2.txt",
        BAM = "{OUTDIR}/{SAMPLE}/STARsolo/{RECIPE}/Aligned.sortedByCoord.out.bam"
    output:
        BAM = "{OUTDIR}/{SAMPLE}/STARsolo/{RECIPE}/Aligned.sortedByCoord.dedup.out.bam"
    threads:
        config['CORES']
    log:
        log = "{OUTDIR}/{SAMPLE}/STARsolo/{RECIPE}/dedup.log"
    run:
        tmp_recipe = RECIPE_DICT[wildcards.SAMPLE]

        shell(
            f"""
            bash scripts/bash/split_dedup.sh \
                {input.BAM} \
                {input.BB_WHITELIST} \
                {threads} \
                {output.BAM} \
                $(dirname {output.BAM})/tmp/dedup \
            | tee {log.log}
            """
        )

# Index the deduplicated .bam file
rule umitools_indexDedupBAM:
    input:
        BAM = '{OUTDIR}/{SAMPLE}/STARsolo/{RECIPE}/Aligned.sortedByCoord.dedup.out.bam'
    output:
        BAI = '{OUTDIR}/{SAMPLE}/STARsolo/{RECIPE}/Aligned.sortedByCoord.dedup.out.bam.bai'
    threads:
        config['CORES']
    run:
        shell(
            f"""
            {EXEC['SAMTOOLS']} index -@ {threads} {input.BAM}
            """
        )

# Split .bam file by strand for IGV browsing
rule strand_split_dedup_bam:
    input:
        DEDUPBAM = '{OUTDIR}/{SAMPLE}/STARsolo/{RECIPE}/Aligned.sortedByCoord.dedup.out.bam'
    output:
        FWDBAM = '{OUTDIR}/{SAMPLE}/STARsolo/{RECIPE}/Aligned.sortedByCoord.dedup.out.fwd.bam',
        REVBAM = '{OUTDIR}/{SAMPLE}/STARsolo/{RECIPE}/Aligned.sortedByCoord.dedup.out.rev.bam'
    threads:
        1
    run:
        shell(
            f"""
            {EXEC['SAMTOOLS']} view -b -F 0x10 {input.DEDUPBAM} > {output.FWDBAM}
            {EXEC['SAMTOOLS']} view -b -f 0x10 {input.DEDUPBAM} > {output.REVBAM}
            """
        )

# Index the split/deduped bam files
rule indexSplitBAMs:
    input:
        FWDBAM = '{OUTDIR}/{SAMPLE}/STARsolo/{RECIPE}/Aligned.sortedByCoord.dedup.out.fwd.bam',
        REVBAM = '{OUTDIR}/{SAMPLE}/STARsolo/{RECIPE}/Aligned.sortedByCoord.dedup.out.rev.bam'
    output:
        FWDBAI = '{OUTDIR}/{SAMPLE}/STARsolo/{RECIPE}/Aligned.sortedByCoord.dedup.out.fwd.bam.bai',
        REVBAI = '{OUTDIR}/{SAMPLE}/STARsolo/{RECIPE}/Aligned.sortedByCoord.dedup.out.rev.bam.bai'
    threads:
        config['CORES']
    run:
        shell(
            f"""
            {EXEC['SAMTOOLS']} index -@ {threads} {input.FWDBAM}
            {EXEC['SAMTOOLS']} index -@ {threads} {input.REVBAM}
            """
        )