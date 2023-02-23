# Remove reads that don't have a corrected spot/cell barcode with samtools, then remove duplicates w/ **umi-tools**
## High mem usage? Check here! https://umi-tools.readthedocs.io/en/latest/faq.html
#TODO: split bam by strand and by chromosome, then dedup each chr!
#TODO: add exec paths for samtools, bamtools
#TODO: chr-split deduping
rule umitools_dedupBAM:
    input:
        BB_WHITELIST = "{OUTDIR}/{sample}/bb/whitelist.txt",
        BB_1 = "{OUTDIR}/{sample}/bb/whitelist_1.txt",
        BB_2 = "{OUTDIR}/{sample}/bb/whitelist_2.txt",
        SORTEDBAM = '{OUTDIR}/{sample}/STARsolo/Aligned.sortedByCoord.out.bam'
    output:
        DEDUPBAM = '{OUTDIR}/{sample}/STARsolo/Aligned.sortedByCoord.dedup.out.bam'
        # TMPBAM = temp('{OUTDIR}/{sample}/tmp/Aligned.sortedByCoord.CBfiltered.out.bam')
    params:
        # OUTPUT_PREFIX='{OUTDIR}/{sample}/umitools_dedup/{sample}'
        # TMPBAM = '{OUTDIR}/{sample}/tmp.bam'
    threads:
        config['CORES']
        #1
    log:
        '{OUTDIR}/{sample}/dedup.log'
    run:
        tmp_chemistry = CHEM_DICT[wildcards.sample]

        #param handling for different alignment strategies
        if "noTrim" in tmp_chemistry:
            whitelist = f"{input.BB_1} {input.BB_2}" #TODO: pretty sure this won't work..
        elif "internalTrim" in tmp_chemistry:
            whitelist = input.BB_WHITELIST
        else:
            whitelist = input.BB_WHITELIST

        # shell(
        #     f"""
        #     {SAMTOOLS_EXEC} view -1 -b \
        #     -@ {threads} \
        #     --tag-file CB:{whitelist} \
        #     -F UB:Z:- \
        #     {input.SORTEDBAM} \
        #     > {output.TMPBAM}

        #     {SAMTOOLS_EXEC} index \
        #     -@ {threads} \
        #     {output.TMPBAM}

        #     {UMITOOLS_EXEC} dedup \
        #     -I {output.TMPBAM} \
        #     --extract-umi-method=tag \
        #     --umi-tag=UB \
        #     --cell-tag=CB \
        #     --method=unique \
        #     --per-cell \
        #     --unmapped-reads=discard \
        #     --output-stats={params.OUTPUT_PREFIX} \
        #     --log {log} \
        #     -S {output.DEDUPBAM}
        #     """
        # )
        shell(
            f"""
            bash scripts/split_dedup.sh {input.SORTEDBAM} {whitelist} {threads} {output.DEDUPBAM} {OUTDIR}/{wildcards.sample}/tmp/dedup | tee {log}
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
    shell:
        """
        {SAMTOOLS_EXEC} index -@ {threads} {input.SORTEDBAM}
        """

#TODO
# rule strand_split_dedup_bam:
#     input:
#         SORTEDBAM = '{OUTDIR}/{sample}/Aligned.sortedByCoord.dedup.out.bam'
#     output:
#         BAI = '{OUTDIR}/{sample}/Aligned.sortedByCoord.dedup.out.bam.bai'
#     threads:
#         config['CORES']
#     shell:
#         """
#         {SAMTOOLS_EXEC} index -@ {threads} {input.SORTEDBAM}
#         """