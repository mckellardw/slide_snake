# Remove reads that don't have a corrected spot/cell barcode with samtools, then remove duplicates w/ **umi-tools**
## High mem usage? Check here! https://umi-tools.readthedocs.io/en/latest/faq.html


# Rule to split BAM files by strand for IGV browsing
rule ilmn_3c_strand_split_bam:
    input:
        BAM="{OUTDIR}/{SAMPLE}/short_read/STARsolo/{RECIPE}/Aligned.sortedByCoord.out.bam",
        BAI="{OUTDIR}/{SAMPLE}/short_read/STARsolo/{RECIPE}/Aligned.sortedByCoord.out.bam.bai",
    output:
        FWDBAM="{OUTDIR}/{SAMPLE}/short_read/STARsolo/{RECIPE}/Aligned.sortedByCoord.out.fwd.bam",
        REVBAM="{OUTDIR}/{SAMPLE}/short_read/STARsolo/{RECIPE}/Aligned.sortedByCoord.out.rev.bam",
    threads: 1
    shell:
        """
        samtools view -b -F 0x10 {input.BAM} > {output.FWDBAM}
        samtools view -b -f 0x10 {input.BAM} > {output.REVBAM}
        """


# Rule to deduplicate forward strand BAM files using bam_dedupByChr.sh
rule ilmn_3c_umitools_dedup_fwdBAM:
    input:
        BAM="{OUTDIR}/{SAMPLE}/short_read/STARsolo/{RECIPE}/Aligned.sortedByCoord.out.fwd.bam",
        BAI="{OUTDIR}/{SAMPLE}/short_read/STARsolo/{RECIPE}/Aligned.sortedByCoord.out.fwd.bam.bai",
    output:
        BAM="{OUTDIR}/{SAMPLE}/short_read/STARsolo/{RECIPE}/Aligned.sortedByCoord.out.fwd.dedup.bam",
    params:
        CELL_TAG="CB",  # corrected cell barcode tag
        UMI_TAG="UB",  # corrected UMI tag
    threads: config["CORES"]
    log:
        log="{OUTDIR}/{SAMPLE}/short_read/logs/{RECIPE}/dedup_fwd.log",
        err="{OUTDIR}/{SAMPLE}/short_read/logs/{RECIPE}/dedup_fwd.err",
    conda:
        f"{workflow.basedir}/envs/umi_tools.yml"
    shell:
        """
        bash scripts/bash/bam_dedupByChr.sh \
            --bam {input.BAM} \
            --cores {threads} \
            --outdir $(dirname {output.BAM}) \
            --tmpdir $(dirname {output.BAM})/tmp/dedup_fwd \
            --celltag {params.CELL_TAG} \
            --umitag {params.UMI_TAG} \
        1> {log.log} \
        2> {log.err}
        """


# Rule to deduplicate reverse strand BAM files using bam_dedupByChr.sh
rule ilmn_3c_umitools_dedup_revBAM:
    input:
        BAM="{OUTDIR}/{SAMPLE}/short_read/STARsolo/{RECIPE}/Aligned.sortedByCoord.out.rev.bam",
        BAI="{OUTDIR}/{SAMPLE}/short_read/STARsolo/{RECIPE}/Aligned.sortedByCoord.out.rev.bam.bai",
    output:
        BAM="{OUTDIR}/{SAMPLE}/short_read/STARsolo/{RECIPE}/Aligned.sortedByCoord.out.rev.dedup.bam",
    params:
        CELL_TAG="CB",  # corrected cell barcode tag
        UMI_TAG="UB",  # corrected UMI tag
    threads: config["CORES"]
    log:
        log="{OUTDIR}/{SAMPLE}/short_read/logs/{RECIPE}/dedup_rev.log",
        err="{OUTDIR}/{SAMPLE}/short_read/logs/{RECIPE}/dedup_rev.err",
    conda:
        f"{workflow.basedir}/envs/umi_tools.yml"
    shell:
        """
        bash scripts/bash/bam_dedupByChr.sh \
            --bam {input.BAM} \
            --cores {threads} \
            --outdir $(dirname {output.BAM}) \
            --tmpdir $(dirname {output.BAM})/tmp/dedup_rev \
            --celltag {params.CELL_TAG} \
            --umitag {params.UMI_TAG} \
        1> {log.log} \
        2> {log.err}
        """


# Rule to merge and sort deduplicated BAM files
rule ilmn_3c_merge_dedup_bam:
    input:
        FWDBAM="{OUTDIR}/{SAMPLE}/short_read/STARsolo/{RECIPE}/Aligned.sortedByCoord.out.fwd.dedup.bam",
        REVBAM="{OUTDIR}/{SAMPLE}/short_read/STARsolo/{RECIPE}/Aligned.sortedByCoord.out.rev.dedup.bam",
    output:
        BAM="{OUTDIR}/{SAMPLE}/short_read/STARsolo/{RECIPE}/Aligned.sortedByCoord.out.dedup.bam",
    threads: 1
    shell:
        """
        samtools merge -f -u {output.BAM} {input.FWDBAM} {input.REVBAM} 
        """
        # \
        # | samtools sort -o {output.BAM}
