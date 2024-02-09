#############################################
## Unmapped read analyses
#############################################

# Run fastqc on unmapped reads; switch names because of STAR weirdness
rule fastqc_unmapped:
    input:
        UNMAPPED1 = '{OUTDIR}/{SAMPLE}/STARsolo/{RECIPE}/Unmapped.out.mate1',
        UNMAPPED2 = '{OUTDIR}/{SAMPLE}/STARsolo/{RECIPE}/Unmapped.out.mate2',
    output:
        UNMAPPED1_FQ = '{OUTDIR}/{SAMPLE}/STARsolo/{RECIPE}/Unmapped.out.mate1.fastq.gz',
        UNMAPPED2_FQ = '{OUTDIR}/{SAMPLE}/STARsolo/{RECIPE}/Unmapped.out.mate2.fastq.gz',
        FQC_DIR = directory('{OUTDIR}/{SAMPLE}/fastqc/unmapped/{RECIPE}')
    params:
        FASTQC_ADAPTERS = config['FASTQC_ADAPTERS']
    threads:
        config['CORES']
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
        UNMAPPED2_FQ = '{OUTDIR}/{SAMPLE}/STARsolo/Unmapped.out.mate2.fastq.gz'
    output:
        BLAST_R2 = '{OUTDIR}/{SAMPLE}/Unmapped.out.mate2_blastResults.txt'
    threads:
        config['CORES']
    params:
        BLASTDB = config['BLASTDB'],
        TMP_FA = '{OUTDIR}/{SAMPLE}/Unmapped.out.mate2.fa',
        N_READS = 1000
    run:
        shell(
            f"""
            zcat {input.UNMAPPED2_FQ} | sed -n '1~4s/^@/>/p;2~4p' > {params.TMP_FA}

            echo "Number of unmapped reads: "
            grep -c ">" {params.TMP_FA}

            {EXEC['VSEARCH']} \
                --sortbysize {params.TMP_FA} \
                --topn {params.N_READS} \
                --output {OUTDIR}/{wildcards.SAMPLE}/top_1000.fa

            {EXEC['BLASTN']} -db {params.BLASTDB}/nt \
                -query {OUTDIR}/{wildcards.SAMPLE}/top_1000.fa \
                -out {output.BLAST_R2} \
                -outfmt '6 qseqid sseqid stitle pident length mismatch gapopen qstart qend sstart send evalue bitscore' \
                -max_target_seqs 5 \
                -num_threads {threads}

            rm {params.TMP_FA}
    		"""
        )
# mv {OUTDIR}/{wildcards.SAMPLE}/Unmapped.out.mate2_blastResults.txt {OUTDIR}/{wildcards.SAMPLE}/Unmapped.out.mate2_blastResults.tsv

# cat {input.UNMAPPED1_FQ} | awk '{{if(NR%4==1) {{printf(">%s\n",substr($0,2));}} else if(NR%4==2) print;}}' > {params.TMP_FA}
