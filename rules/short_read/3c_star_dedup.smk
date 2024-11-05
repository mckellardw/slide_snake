# Remove reads that don't have a corrected spot/cell barcode with samtools, then remove duplicates w/ **umi-tools**
## High mem usage? Check here! https://umi-tools.readthedocs.io/en/latest/faq.html
# TODO: dedup after strand-splitting?
rule ilmn_3c_umitools_dedupBAM:
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
rule ilmn_3c_strand_split_dedup_bam:
    input:
        BAM="{OUTDIR}/{SAMPLE}/STARsolo/short_read/{RECIPE}/Aligned.sortedByCoord.dedup.out.bam",
        BAI="{OUTDIR}/{SAMPLE}/STARsolo/short_read/{RECIPE}/Aligned.sortedByCoord.dedup.out.bam.bai",
    output:
        FWDBAM="{OUTDIR}/{SAMPLE}/STARsolo/short_read/{RECIPE}/Aligned.sortedByCoord.dedup.out.fwd.bam",
        REVBAM="{OUTDIR}/{SAMPLE}/STARsolo/short_read/{RECIPE}/Aligned.sortedByCoord.dedup.out.rev.bam",
    threads: 1
    shell:
        """
        samtools view -b -F 0x10 {input.BAM} > {output.FWDBAM}
        samtools view -b -f 0x10 {input.BAM} > {output.REVBAM}
        """
