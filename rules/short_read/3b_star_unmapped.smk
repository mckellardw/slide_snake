#############################################
## Unmapped read analyses
#############################################


# Run fastqc on unmapped reads; switch names because of STAR weirdness
rule ilmn_3b_fastqc_unmapped:
    input:
        UNMAPPED1="{OUTDIR}/{SAMPLE}/short_read/STARsolo/{RECIPE}/Unmapped.out.mate1",
        UNMAPPED2="{OUTDIR}/{SAMPLE}/short_read/STARsolo/{RECIPE}/Unmapped.out.mate2",
    output:
        UNMAPPED1_FQ="{OUTDIR}/{SAMPLE}/short_read/STARsolo/{RECIPE}/Unmapped.out.mate1.fq.gz",
        UNMAPPED2_FQ="{OUTDIR}/{SAMPLE}/short_read/STARsolo/{RECIPE}/Unmapped.out.mate2.fq.gz",
        fastqcDir=directory("{OUTDIR}/{SAMPLE}/short_read/fastqc/unmapped/{RECIPE}"),
    params:
        FASTQC_ADAPTERS=config["FASTQC_ADAPTERS"],
    # resources:
    threads: config["CORES"]
    conda:
        f"{workflow.basedir}/envs/fastqc.yml"
    shell:
        """        
        mv {input.UNMAPPED1} {input.UNMAPPED2}.fq
        mv {input.UNMAPPED2} {input.UNMAPPED1}.fq

        pigz -p{threads} -f {input.UNMAPPED1}.fq {input.UNMAPPED2}.fq

        mkdir -p {output.fastqcDir}

        fastqc \
            -o {output.fastqcDir} \
            -t {threads} \
            -a {params.FASTQC_ADAPTERS} \
            {output.UNMAPPED1_FQ} {output.UNMAPPED2_FQ}
        """


# Only BLAST R2, which contains the insert (converts .fq to .fa, then removes the .fa file)
rule ilmn_3b_blast_unmapped:
    input:
        UNMAPPED2_FQ="{OUTDIR}/{SAMPLE}/short_read/STARsolo/{RECIPE}/Unmapped.out.mate2.fq.gz",
    output:
        BLAST_R2="{OUTDIR}/{SAMPLE}/unmapped/{RECIPE}/blast/Unmapped.out.mate2_blastResults.txt",
        TMP_FA=temp("{OUTDIR}/{SAMPLE}/unmapped/{RECIPE}/blast/Unmapped.out.mate2.fa"),
        TOP_FA="{OUTDIR}/{SAMPLE}/unmapped/{RECIPE}/blast/top_seqs.fa",
    params:
        BLASTDB=config["BLASTDB"],
        TOPN=10000,
        OUTFMT="6 qseqid sseqid stitle pident length mismatch gapopen qstart qend sstart send evalue bitscore",
    log:
        log="{OUTDIR}/{SAMPLE}/unmapped/{RECIPE}/blast/blast_unmapped.log",
    # resources:
    threads: config["CORES"]
    conda:
        f"{workflow.basedir}/envs/blast.yml"
    shell:
        """
        mkdir -p $(dirname {output.BLAST_R2})

        zcat {input.UNMAPPED2_FQ} \
        | sed -n '1~4s/^@/>/p;2~4p' \
        > {output.TMP_FA}

        echo "Number of unmapped reads: " >> {log.log}
        grep -c ">" {output.TMP_FA}       >> {log.log}

        vsearch \
            --sortbysize {output.TMP_FA} \
            --topn {params.TOPN} \
            --output {output.TOP_FA}

        blastn \
            -db {params.BLASTDB}/nt \
            -query {output.TOP_FA} \
            -out {output.BLAST_R2} \
            -outfmt '{params.OUTFMT}' \
            -max_target_seqs 5 \
            -num_threads {threads}
        """


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
#         UNMAPPED2_FQ = '{OUTDIR}/{SAMPLE}/short_read/STARsolo/{RECIPE}/Unmapped.out.mate2.fq.gz'
#     output:
#         BAM1 = temp('{OUTDIR}/{SAMPLE}/unmapped/{RECIPE}/aligned.bam'), #temp()?
#         BAM2 = '{OUTDIR}/{SAMPLE}unmapped/{RECIPE}/aligned_sorted.bam', #temp()?
#         R2_FQ_BWA_FILTERED  = '{OUTDIR}/{SAMPLE}/unmapped/{RECIPE}/noRibo_R2.fq'
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
#             bwa-mem2 mem \
#                 -t {threads} \
#                 {BWA_REF} \
#                 {input.R2_FQ} \
#             1> {output.BAM1} \
#             2> {log.log} \
#             samtools sort \
#                 -@ {threads} \
#                 -O BAM \
#                 {output.BAM1} \
#             > {output.BAM2}
#             samtools fastq \
#                 -f 4 \
#                 {output.BAM2} \
#             > {output.R2_FQ_BWA_FILTERED}
#             """
#         )
