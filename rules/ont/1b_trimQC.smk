#############################################
## Read trimming & QC rules
#############################################
# fastqc before trimming
rule ont_fastQC_preTrim:
    input:
        MERGED_FQ = '{OUTDIR}/{SAMPLE}/ont/merged.fq.gz'
    output:
        DIR = directory('{OUTDIR}/{SAMPLE}/fastqc/ont_preAdapterScan')
    params:
        adapters = config['FASTQC_ADAPTERS']
    threads:
        config['CORES']
        # min([config['CORES'],8]) # 8 core max
    run:
        shell(
            f"""
            mkdir -p {output.DIR}

            {EXEC['FASTQC']} \
                --outdir {output.DIR} \
                --threads {threads} \
                -a {params.adapters} \
                {input.MERGED_FQ}
            """
        )

# https://github.com/yfukasawa/LongQC
# rule ont_longQC_preTrim:
#     input:
#         MERGED_FQ = '{OUTDIR}/{SAMPLE}/ont/merged.fq.gz'
#     output:
#         DIR = directory('{OUTDIR}/{SAMPLE}/longqc/ont_preAdapterScan')
#     params:
#         adapters = config['FASTQC_ADAPTERS']
#     threads:
#         config['CORES']
#         # min([config['CORES'],8]) # 8 core max
#     run:
#         shell(
#             f"""
#             mkdir -p {output.DIR}

#             {EXEC['LONGQC']} \
#                 -o {output.DIR} \
#                 -p {threads} \
#                 -x ont-ligation \
#                 {input.MERGED_FQ}
#             """
#         )

rule ont_cutadapt:
    input:
        R1_FQ = "{OUTDIR}/{SAMPLE}/ont/adapter_scan_readids/full_len_R1.fq.gz",
        # R1_FQ_Trimmed = "{OUTDIR}/{SAMPLE}/ont/adapter_scan_readids/full_len_internalTrim_R1.fq.gz",
        R2_FQ = "{OUTDIR}/{SAMPLE}/ont/adapter_scan_readids/full_len_R2.fq.gz"
    output:
        R1_FQ = temp('{OUTDIR}/{SAMPLE}/ont/tmp/cut_R1.fq.gz'),
        R2_FQ = temp('{OUTDIR}/{SAMPLE}/ont/tmp/cut_R2.fq.gz'),
        JSON = '{OUTDIR}/{SAMPLE}/ont/cutadapt.json'
    params:
        # R1_LENGTH = 50,
        ADAPTER_COUNT=4, # number of adapters that can be trimmed from each read
        QUALITY_MIN=5, # Low Q score for ONT...
        MIN_R2_LENGTH = 12,
        OVERLAP = 5,
        HOMOPOLYMER_ERROR_RATE = 0.2, # default error rate is 0.1
        POLYA = "A"*100,
        # THREE_PRIME_R2_POLYG = "G"*100,
        POLYT = "T"*100,
        THREE_PRIME_R2_TSO = "AAGCTGGTATCAACGCAGAGTGAATGGG", # SlideSeq TSO - remove any polyadenylated TSOs
        THREE_PRIME_R2_TXG_TSO = "AAGCAGTGGTATCAACGCAGAGTACATGGG", # 10x TSO - remove any polyadenylated TSOs
        FIVE_PRIME_R2_TSO = "CCCATTCACTCTGCGTTGATACCAGCTT", # rev comp of SlideSeq TSO
        FIVE_PRIME_R2_TXG_TSO = "CCCATGTACTCTGCGTTGATACCACTGCTT", # rev-comp of 10x TSO sequence
        THREE_PRIME_R2_SEEKER_BB_ADAPTER = "TCTTCAGCGTTCCCGAGA", # Adapter between BB1 & BB2 in R1 
        FIVE_PRIME_R2_SEEKER_BB_ADAPTER = "AGAGCCCTTGCGACTTCT" # Reverse of the adapter between BB1 & BB2 in R1 
    threads:
        # min([config['CORES'],8]) # 8 core max
        config['CORES']
    log:
        log = '{OUTDIR}/{SAMPLE}/ont/cutadapt.log'
    run:
        # R1_LENGTH = RECIPE_SHEET["R1.finalLength"][recipe]
        R1_LENGTH=10 #TODO- fix this

        R1 = input.R1_FQ
        # R1 = input.R1_FQ_Trimmed
        R2 = input.R2_FQ

        shell(
            f"""
            mkdir -p $(dirname  {output.R1_FQ})
            
            {EXEC['CUTADAPT']} \
                --minimum-length {R1_LENGTH}:{params.MIN_R2_LENGTH} \
                --quality-cutoff {params.QUALITY_MIN} \
                --overlap {params.OVERLAP} \
                --match-read-wildcards \
                --nextseq-trim=20 \
                --times {params.ADAPTER_COUNT} \
                -A "{params.POLYA};max_error_rate={params.HOMOPOLYMER_ERROR_RATE}" \
                -A "{params.POLYT};max_error_rate={params.HOMOPOLYMER_ERROR_RATE}" \
                -A {params.THREE_PRIME_R2_TSO} \
                -A {params.THREE_PRIME_R2_TXG_TSO} \
                -A {params.THREE_PRIME_R2_SEEKER_BB_ADAPTER} \
                -G {params.FIVE_PRIME_R2_TSO} \
                -G {params.FIVE_PRIME_R2_TXG_TSO} \
                -G {params.FIVE_PRIME_R2_SEEKER_BB_ADAPTER} \
                --pair-filter=any \
                -o {output.R1_FQ} \
                -p {output.R2_FQ} \
                --cores {threads} \
                --json {output.JSON} \
                {R1} {R2} \
            1> {log.log}
            """
        )