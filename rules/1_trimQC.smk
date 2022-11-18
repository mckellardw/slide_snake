#############################################
## Read trimming & QC rules
#############################################

# fastqc before trimming
rule preTrim_FastQC_R2:
    input:
        MERGED_R2_FQ = '{OUTDIR}/{sample}/tmp/{sample}_R2.fq.gz'
    output:
        fastqcDir = directory('{OUTDIR}/{sample}/preTrim_fastqc_R2_out')
    params:
        adapters = config['FASTQC_ADAPTERS']
    threads:
        # config['CORES']
        min([config['CORES'],8]) # 8 core max
    shell:
        """
        mkdir -p {output.fastqcDir}
        cd {output.fastqcDir}

        fastqc \
        --outdir {output.fastqcDir} \
        --threads {threads} \
        -a {params.adapters} \
        {input.MERGED_R2_FQ}
        """

# TSO & polyA trimming
rule cutadapt_R2:
    input:
        TRIMMED_R1_FQ = '{OUTDIR}/{sample}/tmp/{sample}_R1.fq.gz',
        TRIMMED_R2_FQ = '{OUTDIR}/{sample}/tmp/{sample}_R2.fq.gz'
    output:
        FINAL_R1_FQ = temp('{OUTDIR}/{sample}/tmp/{sample}_R1_adapterTrim.fq.gz'),
        FINAL_R2_FQ = temp('{OUTDIR}/{sample}/tmp/{sample}_R2_final.fq.gz')
    params:
        R1_SIZE = 50,
        MIN_R2_SIZE = 16,
        CUTADAPT_EXEC = CUTADAPT_EXEC,
        THREE_PRIME_R2_POLYA = "A"*100,
        THREE_PRIME_R2_POLYG = "G"*100,
        THREE_PRIME_R2_POLYT = "T"*100,
        THREE_PRIME_R2_NEXTERA = "CTGTCTCTTATA", # Nextera sequence
        THREE_PRIME_R2_rcNEXTERA = "TATAAGAGACAG", # Rev Comp of Nextera sequence
        THREE_PRIME_R2_TSO = "AAGCTGGTATCAACGCAGAGTGAATGGG", # SlideSeq TSO - remove any polyadenylated TSOs
        THREE_PRIME_R2_ILLUMINA_UNI = "AGATCGGAAGAG", # Illumina Universal
        FIVE_PRIME_R2_TSO = "CCCATTCACTCTGCGTTGATACCAGCTT" # rev comp of SlideSeq TSO
    threads:
        config['CORES']
        # min([config['CORES'],8]) # 8 core max based on recommendations from trim_galore authors
    log:
        log = '{OUTDIR}/{sample}/cutadapt.log'
    shell:
        """
        {params.CUTADAPT_EXEC} \
        --minimum-length {params.R1_SIZE}:{params.MIN_R2_SIZE} \
        --quality-cutoff 20 \
        --overlap 3 \
        --match-read-wildcards \
        --nextseq-trim=20 \
        -A {params.THREE_PRIME_R2_POLYA} \
        -A {params.THREE_PRIME_R2_POLYT} \
        -A {params.THREE_PRIME_R2_TSO} \
        -A {params.THREE_PRIME_R2_NEXTERA} \
        -A {params.THREE_PRIME_R2_rcNEXTERA} \
        -A {params.THREE_PRIME_R2_ILLUMINA_UNI} \
        -G {params.FIVE_PRIME_R2_TSO} \
        --pair-filter=any \
 		-o {output.FINAL_R1_FQ} \
        -p {output.FINAL_R2_FQ} \
        --cores {threads} \
        {input.TRIMMED_R1_FQ} {input.TRIMMED_R2_FQ} 1> {log.log}
        """
        # -A {params.THREE_PRIME_R2_POLYG}X \

# R1 trimming to remove the linker sequence
## Source: https://unix.stackexchange.com/questions/510164/remove-and-add-sequence-information-at-specific-position-in-a-file
rule removeLinker_R1:
    input:
        R1_FQ = '{OUTDIR}/{sample}/tmp/{sample}_R1_adapterTrim.fq.gz'
    output:
        FINAL_R1_FQ = temp('{OUTDIR}/{sample}/tmp/{sample}_R1_final.fq.gz')
    params:
        script = "scripts/linkerRemove_R1.awk",
        CB1end = 8,
        CB2start = 27,
        CB2end = 42
    threads:
        config['CORES']
    run:
        shell(
            f"""
            zcat {input.R1_FQ} | \
            awk -v s={params.CB1end} -v S={params.CB2start} -v E={params.CB2end} -f {params.script} > {OUTDIR}/{wildcards.sample}/tmp/{wildcards.sample}_R1_final.fq

            pigz -p{threads} {OUTDIR}/{wildcards.sample}/tmp/{wildcards.sample}_R1_final.fq
            """
        )

# fastqc on R1 after linker removal & R2 trimming/filtering
rule postTrim_FastQC_R1:
    input:
        FINAL_R1_FQ =  '{OUTDIR}/{sample}/tmp/{sample}_R1_final.fq.gz'
    output:
        fastqcDir = directory('{OUTDIR}/{sample}/postTrim_fastqc_R1_out'),
        # fastqcReport = ''
    threads:
        # config['CORES']
        min([config['CORES'],8]) # 8 core max
    params:
        adapters = config['FASTQC_ADAPTERS']
    shell:
        """
        mkdir -p {output.fastqcDir}

        fastqc \
        --outdir {output.fastqcDir} \
        --threads {threads} \
        -a {params.adapters} \
        {input.FINAL_R1_FQ}
        """

# fastqc after trimming on R2
rule postTrim_FastQC_R2:
    input:
        FINAL_R2_FQ =  '{OUTDIR}/{sample}/tmp/{sample}_R2_final.fq.gz'
    output:
        fastqcDir = directory('{OUTDIR}/{sample}/postTrim_fastqc_R2_out'),
        # fastqcReport = ''
    threads:
        # config['CORES']
        min([config['CORES'],8]) # 8 core max
    params:
        adapters = config['FASTQC_ADAPTERS']
    shell:
        """
        mkdir -p {output.fastqcDir}

        fastqc \
        --outdir {output.fastqcDir} \
        --threads {threads} \
        -a {params.adapters} \
        {input.FINAL_R2_FQ}
        """
