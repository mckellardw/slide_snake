# Remove reads that don't have a corrected spot/cell barcode with samtools, then remove duplicates w/ **umi-tools**
## High mem usage? Check here! https://umi-tools.readthedocs.io/en/latest/faq.html
# TODO: dedup after strand-splitting?
# TODO: add exec paths for samtools, bamtools
rule umitools_dedupBAM:
    input:
        WHITELIST=lambda w: get_whitelist(w, return_type="list"),
        BAM="{OUTDIR}/{SAMPLE}/STARsolo/short_read/{RECIPE}/Aligned.sortedByCoord.out.bam",
    output:
        BAM="{OUTDIR}/{SAMPLE}/STARsolo/short_read/{RECIPE}/Aligned.sortedByCoord.dedup.out.bam",
    params:
        WHITELIST=lambda w: get_whitelist(w),
    threads: config["CORES"]
    log:
        log="{OUTDIR}/{SAMPLE}/STARsolo/short_read/{RECIPE}/dedup.log",
    conda:
        f"{workflow.basedir}/envs/umi_tools.yml"
    shell:
        """
        bash scripts/bash/split_dedup.sh \
            {input.BAM} \
            {input.WHITELIST} \
            {threads} \
            {output.BAM} \
            $(dirname {output.BAM})/tmp/dedup  \
        | tee {log.log}
        """


# Split .bam file by strand for IGV browsing
rule strand_split_dedup_bam:
    input:
        BAM="{OUTDIR}/{SAMPLE}/STARsolo/short_read/{RECIPE}/Aligned.sortedByCoord.dedup.out.bam",
        BAI="{OUTDIR}/{SAMPLE}/STARsolo/short_read/{RECIPE}/Aligned.sortedByCoord.dedup.out.bam.bai",
    output:
        FWDBAM="{OUTDIR}/{SAMPLE}/STARsolo/short_read/{RECIPE}/Aligned.sortedByCoord.dedup.out.fwd.bam",
        REVBAM="{OUTDIR}/{SAMPLE}/STARsolo/short_read/{RECIPE}/Aligned.sortedByCoord.dedup.out.rev.bam",
    threads: 1
    run:
        shell(
            f"""
            {EXEC['SAMTOOLS']} view -b -F 0x10 {input.BAM} > {output.FWDBAM}
            {EXEC['SAMTOOLS']} view -b -f 0x10 {input.BAM} > {output.REVBAM}
            """
        )


# Index the split/deduped bam files
# rule indexSplitBAMs:
#     input:
#         FWDBAM="{OUTDIR}/{SAMPLE}/STARsolo/short_read/{RECIPE}/Aligned.sortedByCoord.dedup.out.fwd.bam",
#         REVBAM="{OUTDIR}/{SAMPLE}/STARsolo/short_read/{RECIPE}/Aligned.sortedByCoord.dedup.out.rev.bam",
#     output:
#         FWDBAI="{OUTDIR}/{SAMPLE}/STARsolo/short_read/{RECIPE}/Aligned.sortedByCoord.dedup.out.fwd.bam.bai",
#         REVBAI="{OUTDIR}/{SAMPLE}/STARsolo/short_read/{RECIPE}/Aligned.sortedByCoord.dedup.out.rev.bam.bai",
#     threads: config["CORES"]
#     run:
#         shell(
#             f"""
#             {EXEC['SAMTOOLS']} index -@ {threads} {input.FWDBAM}
#             {EXEC['SAMTOOLS']} index -@ {threads} {input.REVBAM}
#             """
#         )
