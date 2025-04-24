# Merge .fastq files (in case more than one sequencing run was performed)
rule ilmn_1a_merge_fastqs:
    output:
        MERGED_R1_FQ=temp("{OUTDIR}/{SAMPLE}/short_read/tmp/merged_R1.fq.gz"),
        MERGED_R2_FQ=temp("{OUTDIR}/{SAMPLE}/short_read/tmp/merged_R2.fq.gz"),
    params:
        R1_LIST=lambda wildcards: " ".join(R1_FQS[wildcards.SAMPLE]),
        R2_LIST=lambda wildcards: " ".join(R2_FQS[wildcards.SAMPLE]),
        TMPDIR="{OUTDIR}/{SAMPLE}/short_read/tmp",
    resources:
        mem="16G",
    threads: 1
    log:
        log="{OUTDIR}/{SAMPLE}/short_read/logs/1a_merge_fastqs.log",
        err="{OUTDIR}/{SAMPLE}/short_read/logs/1a_merge_fastqs.err",
    shell:
        """
        bash scripts/bash/merge_formats_short_read.sh \
            -d "{params.TMPDIR}" \
            -1 "{params.R1_LIST}" \
            -2 "{params.R2_LIST}" \
            -o "{output.MERGED_R1_FQ}" \
            -p "{output.MERGED_R2_FQ}" \
            -t "{threads}" \
        1> {log.log} \
        2> {log.err}
        """
