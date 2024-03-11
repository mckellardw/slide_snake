# convert fastq into an unaligned bam for simple passing to alignment tools
rule fq2bam:
    input:
        FINAL_R1_FQ = '{OUTDIR}/{SAMPLE}/tmp/cut_R1.fq.gz',
        FINAL_R2_FQ = '{OUTDIR}/{SAMPLE}/tmp/cut_R2_.fq.gz'
    output:
        uBAM = temp('{OUTDIR}/{SAMPLE}/tmp/unaligned.bam')
    threads:
        config['CORES']
    run:
        shell(
            f"""
            {EXEC['SAMTOOLS']} import \
                -1 {input.FINAL_R1_FQ} \
                -2 {input.FINAL_R2_FQ} \
                -o {output.uBAM}
            """
        )

#TODO: Rule to build whitelist if there isn't one explicitly given
# rule build_whitelist:
#     output:
#         BB_WHITELIST = "{OUTDIR}/{SAMPLE}/bb/whitelist.txt"
#     threads:
#         config['CORES']
#     params:
#     run:
#         shell(
#             f"""
#             {EXEC['UMITOOLS']} whitelist --stdin hgmm_100_R1.fastq.gz \
#                 --bc-pattern=CCCCCCCCCCCCCCCCNNNNNNNNNN \
#                 --set-cell-number=100 \
#                 --log2stderr > whitelist.txt;
#             """
#         )

# rule barcode_ubam:
#     input:
#         uBAM = '{OUTDIR}/{SAMPLE}/tmp/unaligned.bam'
#     output:
#         uBAM = '{OUTDIR}/{SAMPLE}/tmp/unaligned_bc.bam'
#     threads:
#         config['CORES']
#     params:
#     run:
#         shell(
#             f"""
#             #TODO
#             """
#         )