# Merge .fastq files (in case more than one sesquencing run was performed)
rule ilmn_1a_merge_fastqs:
    input:
        R1_FQ=lambda wildcards: R1_FQS[wildcards.SAMPLE],
        R2_FQ=lambda wildcards: R2_FQS[wildcards.SAMPLE],
    output:
        MERGED_R1_FQ=temp("{OUTDIR}/{SAMPLE}/short_read/tmp/merged_R1.fq.gz"),
        MERGED_R2_FQ=temp("{OUTDIR}/{SAMPLE}/short_read/tmp/merged_R2.fq.gz"),
    resources:
        mem="16G",
    threads: 1
    run:
        if (
            len(input.R1_FQ) == 1 & len(input.R2_FQ) == 1
        ):  # shell for single fastq input
            shell(f"cp {input.R1_FQ[0]} {output.MERGED_R1_FQ}")
            shell(f"cp {input.R2_FQ[0]} {output.MERGED_R2_FQ}")
        else:  # multi-fastq input; concatenate inputs
            print("Concatenating", len(input.R1_FQ), ".fastq's for", wildcards.SAMPLE)
            shell(
                f"""
                zcat {" ".join(input.R1_FQ)} > {output.MERGED_R1_FQ.rstrip('.gz')}
                zcat {" ".join(input.R2_FQ)} > {output.MERGED_R2_FQ.rstrip('.gz')}
                pigz -p {threads} {output.MERGED_R1_FQ.rstrip('.gz')} {output.MERGED_R2_FQ.rstrip('.gz')}
                """
            )
