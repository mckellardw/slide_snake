#############################################
## Unmapped read analyses
#############################################


# Run fastqc on unmapped reads; switch names because of STAR weirdness
rule fastqc_unmapped:
    input:
        UNMAPPED1="{OUTDIR}/{SAMPLE}/STARsolo/short_read/{RECIPE}/Unmapped.out.mate1",
        UNMAPPED2="{OUTDIR}/{SAMPLE}/STARsolo/short_read/{RECIPE}/Unmapped.out.mate2",
    output:
        UNMAPPED1_FQ="{OUTDIR}/{SAMPLE}/STARsolo/short_read/{RECIPE}/Unmapped.out.mate1.fastq.gz",
        UNMAPPED2_FQ="{OUTDIR}/{SAMPLE}/STARsolo/short_read/{RECIPE}/Unmapped.out.mate2.fastq.gz",
        FQC_DIR=directory("{OUTDIR}/{SAMPLE}/fastqc/unmapped/{RECIPE}"),
    params:
        FASTQC_ADAPTERS=config["FASTQC_ADAPTERS"],
    threads: config["CORES"]
    run:
        shell(
            f"""
            mv {input.UNMAPPED1} {input.UNMAPPED2}.fastq
            mv {input.UNMAPPED2} {input.UNMAPPED1}.fastq

            {EXEC['PIGZ']} -p{threads} -f {input.UNMAPPED1}.fastq {input.UNMAPPED2}.fastq

            mkdir -p {output.FQC_DIR}

            {EXEC['FASTQC']} \
                -o {output.FQC_DIR} \
                -t {threads} \
                -a {params.FASTQC_ADAPTERS} \
                {output.UNMAPPED1_FQ} {output.UNMAPPED2_FQ}
            """
        )


# Only BLAST R2, which contains the insert (converts .fq to .fa, then removes the .fa file)
## TODO: change demux step to fastx-collapser
rule blast_unmapped:
    input:
        UNMAPPED2_FQ="{OUTDIR}/{SAMPLE}/STARsolo/short_read/{RECIPE}/Unmapped.out.mate2.fastq.gz",
    output:
        BLAST_R2="{OUTDIR}/{SAMPLE}/unmapped/{RECIPE}/blast/Unmapped.out.mate2_blastResults.txt",
        TMP_FA=temp("{OUTDIR}/{SAMPLE}/unmapped/{RECIPE}/blast/Unmapped.out.mate2.fa"),
        TOP_FA="{OUTDIR}/{SAMPLE}/unmapped/{RECIPE}/blast/top_seqs.fa",
    threads: config["CORES"]
    params:
        BLASTDB=config["BLASTDB"],
        TOPN=10000,
        OUTFMT="6 qseqid sseqid stitle pident length mismatch gapopen qstart qend sstart send evalue bitscore",
    log:
        log="{OUTDIR}/{SAMPLE}/unmapped/{RECIPE}/blast/blast_unmapped.log",
    run:
        shell(
            f"""
            mkdir -p $(dirname {output.BLAST_R2})

            zcat {input.UNMAPPED2_FQ} \
            | sed -n '1~4s/^@/>/p;2~4p' \
            > {output.TMP_FA}

            echo "Number of unmapped reads: " >> {log.log}
            grep -c ">" {output.TMP_FA}       >> {log.log}

            {EXEC['VSEARCH']} \
                --sortbysize {output.TMP_FA} \
                --topn {params.TOPN} \
                --output {output.TOP_FA}

            {EXEC['BLASTN']} \
                -db {params.BLASTDB}/nt \
                -query {output.TOP_FA} \
                -out {output.BLAST_R2} \
                -outfmt '{params.OUTFMT}' \
                -max_target_seqs 5 \
                -num_threads {threads}
            """
        )


# mv {OUTDIR}/{wildcards.SAMPLE}/Unmapped.out.mate2_blastResults.txt {OUTDIR}/{wildcards.SAMPLE}/Unmapped.out.mate2_blastResults.tsv

# cat {input.UNMAPPED1_FQ} | awk '{{if(NR%4==1) {{printf(">%s\n",substr($0,2));}} else if(NR%4==2) print;}}' > {params.TMP_FA}

# rule bwa_index_phix:
#     input:
#        PHIX_IDX = 'resources/phix/bwa.idx',
#     output:
#     params:
#         MEMLIMIT = config['MEMLIMIT']
#     # log:
#     #     log = ""
#     threads:
#         config['CORES']
#     run:
#         shell(
#             f"""

#             """
#         )
# TODO
# temp rule - just checking phiX contamination
# rule unmapped_phix_bwa:
#     input:
#         PHIX_IDX = 'resources/phix/bwa.idx',
#         UNMAPPED2_FQ = '{OUTDIR}/{SAMPLE}/STARsolo/short_read/{RECIPE}/Unmapped.out.mate2.fastq.gz'
#     output:
#         BAM1 = temp('{OUTDIR}/{SAMPLE}/unmapped/{RECIPE}/aligned.bam'), #temp()?
#         BAM2 = '{OUTDIR}/{SAMPLE}unmapped/{RECIPE}/aligned_sorted.bam', #temp()?
#         R2_FQ_BWA_FILTERED  = '{OUTDIR}/{SAMPLE}/unmapped/{RECIPE}/final_filtered_R2.fq'
#     params:
#         MEMLIMIT = config['MEMLIMIT']
#     log:
#         log = "{OUTDIR}/{SAMPLE}/unmapped/{RECIPE}/bwa_phix.log"
#     threads:
#         config['CORES']
#     run:
#         shell(
#             f"""
#             mkdir -p $(dirname {output.BAM1})
#             {EXEC['BWA']} mem \
#                 -t {threads} \
#                 {BWA_REF} \
#                 {input.R2_FQ} \
#             1> {output.BAM1} \
#             2> {log.log} \
#             {EXEC['SAMTOOLS']} sort \
#                 -@ {threads} \
#                 -O BAM \
#                 {output.BAM1} \
#             > {output.BAM2}
#             {EXEC['SAMTOOLS']} fastq \
#                 -f 4 \
#                 {output.BAM2} \
#             > {output.R2_FQ_BWA_FILTERED}
#             """
#         )
