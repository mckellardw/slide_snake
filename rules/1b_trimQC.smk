#############################################
## Read trimming & QC rules
#############################################
# fastqc on R2 before trimming
rule preTrim_FastQC_R1:
    input:
        MERGED_R1_FQ = '{OUTDIR}/{sample}/tmp/{sample}_R1.fq.gz'
    output:
        fastqcDir = directory('{OUTDIR}/{sample}/fastqc_preTrim_R1')
    params:
        adapters = config['FASTQC_ADAPTERS']
    threads:
        config['CORES']
        # min([config['CORES'],8]) # 8 core max
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

# fastqc on R2 before trimming
rule preTrim_FastQC_R2:
    input:
        MERGED_R2_FQ = '{OUTDIR}/{sample}/tmp/{sample}_R2.fq.gz'
    output:
        fastqcDir = directory('{OUTDIR}/{sample}/fastqc_preTrim_R2')
    params:
        ADAPTERS = config['FASTQC_ADAPTERS']
    threads:
        config['CORES']
        # min([config['CORES'],8]) # 8 core max
    run:
        shell(
            f"""
            mkdir -p {output.fastqcDir}

            {FASTQC_EXEC} \
            --outdir {output.fastqcDir} \
            --threads {threads} \
            -a {params.ADAPTERS} \
            {input.MERGED_R2_FQ}
            """
        )


#TODO merge trimming rules
# Internal adapter trimming on R1
# rule internal_adapter_trim_R1:
#     input:
#         MERGED_R1_FQ = '{OUTDIR}/{sample}/tmp/{sample}_R1.fq.gz',
#         # MERGED_R2_FQ = '{OUTDIR}/{sample}/tmp/{sample}_R2.fq.gz'
#     output:
#         INTERNAL_R1_FQ = temp('{OUTDIR}/{sample}/tmp/{sample}_R1_InternalTrim.fq.gz'),
#         INTERNAL_TRIM_QC_LOG = '{OUTDIR}/{sample}/internal_trim_qc.txt'
#         # FINAL_R2_FQ = temp('{OUTDIR}/{sample}/tmp/{sample}_R2_final.fq.gz')
#     params:
#         TMPDIR = "{OUTDIR}/{sample}/tmp/seqkit",
#         INTERNAL_ADAPTER = config["R1_INTERNAL_ADAPTER"] # Curio R1 internal adapter
#     threads:
#         config['CORES']        
#     log:
#         '{OUTDIR}/{sample}/internal_adapter_trim_R1.log'
#     run:
#         shell(
#             f"""
#             python scripts/internal_adapter_trim_R1.py {params.INTERNAL_ADAPTER} {output.INTERNAL_TRIM_QC_LOG} {threads} {params.TMPDIR} {input.MERGED_R1_FQ} {output.INTERNAL_R1_FQ} | tee {log}
#             """
#         )

# R1 trimming to remove the linker sequence
## Source: https://unix.stackexchange.com/questions/510164/remove-and-add-sequence-information-at-specific-position-in-a-file
# rule removeLinker_R1:
#     input:
#         R1_FQ = '{OUTDIR}/{sample}/tmp/{sample}_R1.fq.gz'
#     output:
#         FINAL_R1_FQ = temp('{OUTDIR}/{sample}/tmp/{sample}_R1_HardTrim.fq.gz')
#     params:
#         script = "scripts/hardTrimFq.awk",
#         CB1end = 8, #TODO- move to config? or recipe_sheet?
#         CB2start = 27,
#         CB2end = 42
#     threads:
#         config['CORES']
#     run:
#         # tmp_recipe = RECIPE_DICT[wildcards.sample]
#         shell(
#             f"""
#             zcat {input.R1_FQ} | \
#             awk -v s={params.CB1end} -v S={params.CB2start} -v E={params.CB2end} -f {params.script} > {OUTDIR}/{wildcards.sample}/tmp/{wildcards.sample}_R1_HardTrim.fq

#             pigz -f -p{threads} {OUTDIR}/{wildcards.sample}/tmp/{wildcards.sample}_R1_HardTrim.fq
#             """
#         )

# Trimming for R1 to handle Curio adapter issues. See README for recipe details (#TODO)
rule R1_trimming:
    input:
        R1_FQ = '{OUTDIR}/{sample}/tmp/{sample}_R1.fq.gz'
    output:
        R1_FQ = temp('{OUTDIR}/{sample}/tmp/{sample}_R1_Trimmed.fq.gz')
    params:
        script = "scripts/hardTrimFq.awk",
        CB1end = 8, #TODO- move to config? or recipe_sheet?
        CB2start = 27,
        CB2end = 42,
        INTERNAL_TRIM_QC_LOG = '{OUTDIR}/{sample}/internal_trim_qc.txt',
        TMPDIR = "{OUTDIR}/{sample}/tmp/seqkit",
        INTERNAL_ADAPTER = config["R1_INTERNAL_ADAPTER"] # Curio R1 internal adapter
    threads:
        config['CORES']
    log:
        '{OUTDIR}/{sample}/R1_trimming.log'
    run:
        tmp_recipe = RECIPE_DICT[wildcards.sample]
        R1_LENGTH = RECIPE_SHEET["R1.finalLength"][tmp_recipe]
        R1 = input.R1_FQ

        #param handling for different alignment strategies
        if "hardTrim" in tmp_recipe:
            # R1 = input.R1_FQ_HardTrim
            shell( # "Hard" trimming, to remove the adapter based on hard-coded base positions
                f"""
                zcat {input.R1_FQ} | \
                awk -v s={params.CB1end} -v S={params.CB2start} -v E={params.CB2end} -f {params.script} > {OUTDIR}/{wildcards.sample}/tmp/{wildcards.sample}_R1_Trimmed.fq

                pigz -f -p{threads} {OUTDIR}/{wildcards.sample}/tmp/{wildcards.sample}_R1_Trimmed.fq 

                echo "Hard trimming performed on {R1}" > {log}
                """
            )
        elif "internalTrim" in tmp_recipe:
            #TODO- rewrite/speed up internal trimming!
            # R1 = input.R1_FQ_InternalTrim
            shell( # Internal trimming to cut out the adapter sequence
                f"""
                python scripts/internal_adapter_trim_R1.py {params.INTERNAL_ADAPTER} {params.INTERNAL_TRIM_QC_LOG} {threads} {params.TMPDIR} {R1} {output.R1_FQ} | tee {log}
                """
            )
        else:
            # R1 = input.R1_FQ
            shell( # Rename R1_FQ if no trimming needed
                f"""
                mv {R1} {output.R1_FQ} 
                echo "No trimming performed on {R1}..." {log}
                """
            )

# TSO & polyA trimming
#TODO: add "{ADAPTER};noindels" to adapter sequence trimming? - *Note- do not do this for the BB_ADAPTER
rule cutadapt:
    input:
        # R1_FQ = '{OUTDIR}/{sample}/tmp/{sample}_R1.fq.gz',
        # R1_FQ_HardTrim = '{OUTDIR}/{sample}/tmp/{sample}_R1_HardTrim.fq.gz',
        # R1_FQ_InternalTrim = '{OUTDIR}/{sample}/tmp/{sample}_R1_InternalTrim.fq.gz',
        R1_FQ_Trimmed = '{OUTDIR}/{sample}/tmp/{sample}_R1_Trimmed.fq.gz',
        R2_FQ = '{OUTDIR}/{sample}/tmp/{sample}_R2.fq.gz'
    output:
        FINAL_R1_FQ = temp('{OUTDIR}/{sample}/tmp/{sample}_R1_final.fq.gz'),
        FINAL_R2_FQ = temp('{OUTDIR}/{sample}/tmp/{sample}_R2_final.fq.gz')
    params:
        # R1_LENGTH = 50,
        MIN_R2_LENGTH = 12,
        CUTADAPT_EXEC = CUTADAPT_EXEC,
        OVERLAP = 5,
        HOMOPOLYMER_ERROR_RATE = 0.2, # default error rate is 0.1
        # THREE_PRIME_R1_POLYA = "A"*100,
        THREE_PRIME_R2_POLYA = "A"*100,
        # THREE_PRIME_R2_POLYG = "G"*100,
        THREE_PRIME_R2_POLYT = "T"*100,
        THREE_PRIME_R2_NEXTERA = "CTGTCTCTTATA", # Nextera sequence
        THREE_PRIME_R2_rcNEXTERA = "TATAAGAGACAG", # Rev Comp of Nextera sequence
        THREE_PRIME_R2_TSO = "AAGCTGGTATCAACGCAGAGTGAATGGG", # SlideSeq TSO - remove any polyadenylated TSOs
        THREE_PRIME_R2_TXG_TSO = "AAGCAGTGGTATCAACGCAGAGTACATGGG", # 10x TSO - remove any polyadenylated TSOs
        THREE_PRIME_R2_ILLUMINA_UNI = "AGATCGGAAGAG", # Illumina Universal
        FIVE_PRIME_R2_TSO = "CCCATTCACTCTGCGTTGATACCAGCTT", # rev comp of SlideSeq TSO
        FIVE_PRIME_R2_TXG_TSO = "CCCATGTACTCTGCGTTGATACCACTGCTT", # rev-comp of 10x TSO sequence
        THREE_PRIME_R2_SEEKER_BB_ADAPTER = "TCTTCAGCGTTCCCGAGA", # Adapter between BB1 & BB2 in R1 
        FIVE_PRIME_R2_SEEKER_BB_ADAPTER = "AGAGCCCTTGCGACTTCT" # Reverse of the adapter between BB1 & BB2 in R1 
    threads:
        # min([config['CORES'],8]) # 8 core max
        config['CORES']
    log:
        log = '{OUTDIR}/{sample}/cutadapt.log'
    run:
        tmp_recipe = RECIPE_DICT[wildcards.sample]
        R1_LENGTH = RECIPE_SHEET["R1.finalLength"][tmp_recipe]

        #param handling for different alignment strategies
        # if "noTrim" in tmp_recipe:
        #     R1 = input.R1_FQ
        # elif "internalTrim" in tmp_recipe:
        #     R1 = input.R1_FQ_InternalTrim
        # else:
        #     R1 = input.R1_FQ_HardTrim

        R1 = input.R1_FQ_Trimmed
        R2 = input.R2_FQ

        shell(
            f"""
            {params.CUTADAPT_EXEC} \
            --minimum-length {R1_LENGTH}:{params.MIN_R2_LENGTH} \
            --quality-cutoff 20 \
            --overlap {params.OVERLAP} \
            --match-read-wildcards \
            --nextseq-trim=20 \
            -A "{params.THREE_PRIME_R2_POLYA};max_error_rate={params.HOMOPOLYMER_ERROR_RATE}" \
            -A "{params.THREE_PRIME_R2_POLYT};max_error_rate={params.HOMOPOLYMER_ERROR_RATE}" \
            -A {params.THREE_PRIME_R2_TSO} \
            -A {params.THREE_PRIME_R2_TXG_TSO} \
            -A {params.THREE_PRIME_R2_SEEKER_BB_ADAPTER} \
            -A {params.THREE_PRIME_R2_NEXTERA} \
            -A {params.THREE_PRIME_R2_rcNEXTERA} \
            -A {params.THREE_PRIME_R2_ILLUMINA_UNI} \
            -G {params.FIVE_PRIME_R2_TSO} \
            -G {params.FIVE_PRIME_R2_TXG_TSO} \
            -G {params.FIVE_PRIME_R2_SEEKER_BB_ADAPTER} \
            --pair-filter=any \
     		-o {output.FINAL_R1_FQ} \
            -p {output.FINAL_R2_FQ} \
            --cores {threads} \
            {R1} {R2} 1> {log.log}
            """
        )


# fastqc on R1 after linker removal & R2 trimming/filtering
rule postTrim_FastQC_R1:
    input:
        FINAL_R1_FQ =  '{OUTDIR}/{sample}/tmp/{sample}_R1_final.fq.gz'
    output:
        fastqcDir = directory('{OUTDIR}/{sample}/fastqc_postTrim_R1'),
        # fastqcReport = ''
    threads:
        config['CORES']
        # min([config['CORES'],8]) # 8 core max
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
        fastqcDir = directory('{OUTDIR}/{sample}/fastqc_postTrim_R2'),
        # fastqcReport = ''
    threads:
        config['CORES']
        # min([config['CORES'],8]) # 8 core max
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
            {input.FINAL_R2_FQ}
            """
        )
