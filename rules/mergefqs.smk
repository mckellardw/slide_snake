# Merge .fastq files (in case more than one sesquencing run was performed)
rule merge_fastqs:
    output:
        MERGED_R1_FQ = temp('{OUTDIR}/{sample}/tmp/{sample}_R1.fq.gz'),
        MERGED_R2_FQ = temp('{OUTDIR}/{sample}/tmp/{sample}_R2.fq.gz')
    params:
        TMP_DIR = '{OUTDIR}/{sample}/tmp',
        R1_FQ = lambda wildcards: R1_FQS[wildcards.sample],
        R2_FQ = lambda wildcards: R2_FQS[wildcards.sample]
    threads:
        config['CORES']
    # conda:
    #     "STARsolo"
    run:
        if len(params.R1_FQ.split(" "))==1 & len(params.R2_FQ.split(" "))==1: # shell for single fastq input
            shell("cp {params.R1_FQ} {output.MERGED_R1_FQ}")
            shell("cp {params.R2_FQ} {output.MERGED_R2_FQ}")
        else: # shell enablinging multi-fast input; concatenate inputs
            print("Concatenating",len(params.R1_FQ.split(" ")), ".fastq's for", wildcards.sample)
            shell("mkdir -p {params.TMP_DIR}")
            shell("zcat {params.R1_FQ} > {params.TMP_DIR}/{wildcards.sample}_R1.fq")
            shell("zcat {params.R2_FQ} > {params.TMP_DIR}/{wildcards.sample}_R2.fq")
            shell("pigz -p {threads} {params.TMP_DIR}/*.fq")
