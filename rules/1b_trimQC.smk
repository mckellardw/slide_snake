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
        min([config['CORES'],8]) # 8 core max
    run:
        shell(
            f"""
            mkdir -p {output.fastqcDir}

            {FASTQC_EXEC} \
            --outdir {output.fastqcDir} \
            --threads {threads} \
            -a {params.adapters} \
            {input.MERGED_R2_FQ}
            """
        )

rule preTrim_FastQC_R1:
    input:
        MERGED_R1_FQ = '{OUTDIR}/{sample}/tmp/{sample}_R1.fq.gz'
    output:
        fastqcDir = directory('{OUTDIR}/{sample}/preTrim_fastqc_R1_out')
    params:
        adapters = config['FASTQC_ADAPTERS']
    threads:
        min([config['CORES'],8]) # 8 core max
    run:
        shell(
            f"""
            mkdir -p {output.fastqcDir}

            {FASTQC_EXEC} \
            --outdir {output.fastqcDir} \
            --threads {threads} \
            -a {params.adapters} \
            {input.MERGED_R1_FQ}
            """
        )

# TSO & polyA trimming
#TODO: add "{ADAPTER};noindels" to adapter sequence trimming? - *Note- do not do this for the BB_ADAPTER
rule cutadapt_R2:
    input:
        TRIMMED_R1_FQ = '{OUTDIR}/{sample}/tmp/{sample}_R1.fq.gz',
        TRIMMED_R2_FQ = '{OUTDIR}/{sample}/tmp/{sample}_R2.fq.gz'
    output:
        FINAL_R1_FQ = temp('{OUTDIR}/{sample}/tmp/{sample}_R1_adapterTrim.fq.gz'),
        FINAL_R2_FQ = temp('{OUTDIR}/{sample}/tmp/{sample}_R2_final.fq.gz')
    params:
        R1_SIZE = 50,
        MIN_R2_SIZE = 12,
        CUTADAPT_EXEC = CUTADAPT_EXEC,
        HOMOPOLYMER_ERROR_RATE = 0.2, # default error rate is 0.1
        # THREE_PRIME_R1_POLYA = "A"*100,
        THREE_PRIME_R2_POLYA = "A"*100,
        # THREE_PRIME_R2_POLYG = "G"*100,
        THREE_PRIME_R2_POLYT = "T"*100,
        THREE_PRIME_R2_NEXTERA = "CTGTCTCTTATA", # Nextera sequence
        THREE_PRIME_R2_rcNEXTERA = "TATAAGAGACAG", # Rev Comp of Nextera sequence
        THREE_PRIME_R2_TSO = "AAGCTGGTATCAACGCAGAGTGAATGGG", # SlideSeq TSO - remove any polyadenylated TSOs
        THREE_PRIME_R2_ILLUMINA_UNI = "AGATCGGAAGAG", # Illumina Universal
        FIVE_PRIME_R2_TSO = "CCCATTCACTCTGCGTTGATACCAGCTT", # rev comp of SlideSeq TSO
        THREE_PRIME_R2_SEEKER_BB_ADAPTER = "TCTTCAGCGTTCCCGAGA", # Adapter between BB1 & BB2 in R1 
        FIVE_PRIME_R2_SEEKER_BB_ADAPTER = "TCTTCAGCGTTCCCGAGA" # Reverse of the adapter between BB1 & BB2 in R1 
    threads:
        # min([config['CORES'],8]) # 8 core max
        config['CORES']
    log:
        log = '{OUTDIR}/{sample}/cutadapt.log'
    run:
        shell(
            f"""
            {params.CUTADAPT_EXEC} \
            --minimum-length {params.R1_SIZE}:{params.MIN_R2_SIZE} \
            --quality-cutoff 20 \
            --overlap 3 \
            --match-read-wildcards \
            --nextseq-trim=20 \
            -A "{params.THREE_PRIME_R2_POLYA};max_error_rate={params.HOMOPOLYMER_ERROR_RATE}" \
            -A "{params.THREE_PRIME_R2_POLYT};max_error_rate={params.HOMOPOLYMER_ERROR_RATE}" \
            -A {params.THREE_PRIME_R2_TSO} \
            -A {params.THREE_PRIME_R2_SEEKER_BB_ADAPTER} \
            -A {params.THREE_PRIME_R2_NEXTERA} \
            -A {params.THREE_PRIME_R2_rcNEXTERA} \
            -A {params.THREE_PRIME_R2_ILLUMINA_UNI} \
            -G {params.FIVE_PRIME_R2_TSO} \
            -G {params.FIVE_PRIME_R2_SEEKER_BB_ADAPTER} \
            --pair-filter=any \
     		-o {output.FINAL_R1_FQ} \
            -p {output.FINAL_R2_FQ} \
            --cores {threads} \
            {input.TRIMMED_R1_FQ} {input.TRIMMED_R2_FQ} 1> {log.log}
            """
        )

# Internal adapter trimming on R1
rule internal_adapter_trim_R1:
    input:
        MERGED_R1_FQ = '{OUTDIR}/{sample}/tmp/{sample}_R1_adapterTrim.fq.gz',
        # MERGED_R2_FQ = '{OUTDIR}/{sample}/tmp/{sample}_R2.fq.gz'
    output:
        FINAL_R1_FQ = temp('{OUTDIR}/{sample}/tmp/{sample}_R1_finalInternalTrim.fq.gz'),
        # FINAL_R2_FQ = temp('{OUTDIR}/{sample}/tmp/{sample}_R2_final.fq.gz')
    params:
        TMPDIR = "{OUTDIR}/{sample}/tmp/seqtk",
        INTERNAL_ADAPTER = "TCTTCAGCGTTCCCGAGA" # Curio adapter
    threads:
        config['CORES']        
    log:
        '{OUTDIR}/{sample}/internal_adapter_trim_R1.log'
    run:
        shell(
            f"""
            python scripts/internal_adapter_trim_R1.py {params.INTERNAL_ADAPTER} {log} {threads} {params.TMPDIR} {input.MERGED_R1_FQ} {output.FINAL_R1_FQ}
            """
        )


# R1 trimming to remove the linker sequence
## Source: https://unix.stackexchange.com/questions/510164/remove-and-add-sequence-information-at-specific-position-in-a-file
rule removeLinker_R1:
    input:
        R1_FQ = '{OUTDIR}/{sample}/tmp/{sample}_R1_adapterTrim.fq.gz'
    output:
        FINAL_R1_FQ = temp('{OUTDIR}/{sample}/tmp/{sample}_R1_finalHardTrim.fq.gz')
    params:
        script = "scripts/linkerRemove_R1.awk",
        CB1end = 8, #TODO- move to config!
        CB2start = 27,
        CB2end = 42
    threads:
        config['CORES']
    run:
        # tmp_chemistry = CHEM_DICT[wildcards.sample]
        shell(
            f"""
            zcat {input.R1_FQ} | \
            awk -v s={params.CB1end} -v S={params.CB2start} -v E={params.CB2end} -f {params.script} > {OUTDIR}/{wildcards.sample}/tmp/{wildcards.sample}_R1_finalHardTrim.fq

            pigz -p{threads} {OUTDIR}/{wildcards.sample}/tmp/{wildcards.sample}_R1_finalHardTrim.fq
            """
        )


# fastqc on R1 after linker removal & R2 trimming/filtering
rule postTrim_FastQC_R1:
    input:
        FINAL_R1_FQ =  '{OUTDIR}/{sample}/tmp/{sample}_R1_adapterTrim.fq.gz'
    output:
        fastqcDir = directory('{OUTDIR}/{sample}/postTrim_fastqc_R1_out'),
        # fastqcReport = ''
    threads:
        min([config['CORES'],8]) # 8 core max
    params:
        adapters = config['FASTQC_ADAPTERS']
    run:
        shell(
            f"""
            mkdir -p {output.fastqcDir}

            {FASTQC_EXEC} \
            --outdir {output.fastqcDir} \
            --threads {threads} \
            -a {params.adapters} \
            {input.FINAL_R1_FQ}
            """
        )

# fastqc after trimming on R2
rule postTrim_FastQC_R2:
    input:
        FINAL_R2_FQ =  '{OUTDIR}/{sample}/tmp/{sample}_R2_final.fq.gz'
    output:
        fastqcDir = directory('{OUTDIR}/{sample}/postTrim_fastqc_R2_out'),
        # fastqcReport = ''
    threads:
        min([config['CORES'],8]) # 8 core max
    params:
        adapters = config['FASTQC_ADAPTERS']
    run:
        shell(
        f"""
            mkdir -p {output.fastqcDir}

            fastqc \
            --outdir {output.fastqcDir} \
            --threads {threads} \
            -a {params.adapters} \
            {input.FINAL_R2_FQ}
            """
        )
