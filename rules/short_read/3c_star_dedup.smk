# Remove reads that don't have a corrected spot/cell barcode with samtools, then remove duplicates w/ **umi-tools**
## High mem usage? Check here! https://umi-tools.readthedocs.io/en/latest/faq.html
# TODO: dedup after strand-splitting?

# Rule to deduplicate BAM files using bam_dedupByChr.sh
rule ilmn_3c_umitools_dedupBAM:
    input:
        BAM="{OUTDIR}/{SAMPLE}/short_read/STARsolo/{RECIPE}/Aligned.sortedByCoord.out.bam",
    output:
        BAM="{OUTDIR}/{SAMPLE}/short_read/STARsolo/{RECIPE}/Aligned.sortedByCoord.out.dedup.bam",
    params:
        CELL_TAG="CB",  # default cell tag
        UMI_TAG="UB",   # default UMI tag
    threads: config["CORES"]
    log:
        log="{OUTDIR}/{SAMPLE}/short_read/misc_logs/{RECIPE}/dedup.log",
        err="{OUTDIR}/{SAMPLE}/short_read/misc_logs/{RECIPE}/dedup.err",
    conda:
        f"{workflow.basedir}/envs/umi_tools.yml"
    shell:
        """
        bash scripts/bash/bam_dedupByChr.sh \
            --bam {input.BAM} \
            --cores {threads} \
            --outdir $(dirname {output.BAM}) \
            --tmpdir $(dirname {output.BAM})/tmp/dedup \
            --celltag {params.CELL_TAG} \
            --umitag {params.UMI_TAG} \
        1> {log.log} \
        2> {log.err}
        """

# Rule to split deduplicated BAM files by strand for IGV browsing
rule ilmn_3c_strand_split_dedup_bam:
    input:
        BAM="1IR}/{SAMPLE}/short_read/STARsolo/{RECIPE}/Aligned.sortedByCoord.out.dedup.bam",
        BAI="{OUTDIR}/{SAMPLE}/short_read/STARsolo/{RECIPE}/Aligned.sortedByCoord.dedup.out.bam.bai",
    output:
        FWDBAM="{OUTDIR}/{SAMPLE}/short_read/STARsolo/{RECIPE}/Aligned.sortedByCoord.out.dedup.fwd.bam",
        REVBAM="{OUTDIR}/{SAMPLE}/short_read/STARsolo/{RECIPE}/Aligned.sortedByCoord.out.dedup.rev.bam",
    threads: 1
    shell:
        """
        samtools view -b -F 0x10 {input.BAM} > {output.FWDBAM}
        samtools view -b -f 0x10 {input.BAM} > {output.REVBAM}
        """
