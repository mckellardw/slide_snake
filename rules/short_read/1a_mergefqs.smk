# Merge .fastq files (in case more than one sesquencing run was performed)
rule merge_fastqs:
    output:
        MERGED_R1_FQ=temp("{OUTDIR}/{SAMPLE}/tmp/merged_R1.fq.gz"),
        MERGED_R2_FQ=temp("{OUTDIR}/{SAMPLE}/tmp/merged_R2.fq.gz"),
    params:
        TMP_DIR="{OUTDIR}/{SAMPLE}/tmp",
        R1_FQ=lambda wildcards: R1_FQS[wildcards.SAMPLE],
        R2_FQ=lambda wildcards: R2_FQS[wildcards.SAMPLE],
    resources:
        mem="16G",
    threads: 1
    run:
        if (
            len(params.R1_FQ) == 1 & len(params.R2_FQ) == 1
        ):  # shell for single fastq input
            shell(f"cp {params.R1_FQ[0]} {output.MERGED_R1_FQ}")
            shell(f"cp {params.R2_FQ[0]} {output.MERGED_R2_FQ}")
        else:  # multi-fastq input; concatenate inputs
            print("Concatenating", len(params.R1_FQ), ".fastq's for", wildcards.SAMPLE)
            shell(
                f"""
                mkdir -p {params.TMP_DIR}
                zcat {" ".join(params.R1_FQ)} > {params.TMP_DIR}/merged_R1.fq
                zcat {" ".join(params.R2_FQ)} > {params.TMP_DIR}/merged_R2.fq
                pigz -p {threads} {params.TMP_DIR}/*.fq
                """
            )
