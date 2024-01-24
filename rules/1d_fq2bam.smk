# convert fastq into an unaligned bam for simple passing to alignment tools
rule fq2bam:
    input:
        FINAL_R1_FQ = '{OUTDIR}/{sample}/tmp/{sample}_R1_final.fq.gz',
        FINAL_R2_FQ = '{OUTDIR}/{sample}/tmp/{sample}_R2_final.fq.gz'
    output:
        uBAM = temp('{OUTDIR}/{sample}/tmp/unaligned.bam')
    threads:
        config['CORES']
    params:
    run:
        shell(
            f"""
            {EXEC['SAMTOOLS']} import -i \
                -1 {input.FINAL_R1_FQ} \
                -2 {input.FINAL_R2_FQ} \
                -o {output.uBAM}
            """
        )