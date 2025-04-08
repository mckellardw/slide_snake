# Merge .fastq files (in case more than one sequencing run was performed)
rule ilmn_1a_merge_fastqs:
    input:
        R1_FQ=lambda wildcards: R1_FQS[wildcards.SAMPLE],
        R2_FQ=lambda wildcards: R2_FQS[wildcards.SAMPLE],
    output:
        MERGED_R1_FQ=temp("{OUTDIR}/{SAMPLE}/short_read/tmp/merged_R1.fq.gz"),
        MERGED_R2_FQ=temp("{OUTDIR}/{SAMPLE}/short_read/tmp/merged_R2.fq.gz"),
    params:
        TMPDIR="{OUTDIR}/{SAMPLE}/short_read/tmp",
        READS=lambda wildcards: " ".join(
            R1_FQS[wildcards.SAMPLE] + R2_FQS[wildcards.SAMPLE]
        ),
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
            -r "{params.READS}" \
            -o "{output.MERGED_R1_FQ}" \
            -p "{output.MERGED_R2_FQ}" \
            -t "{threads}" \
        1> {log.log} \
        2> {log.err}
        """
