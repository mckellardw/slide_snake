#############################################
## Read trimming & QC rules
#############################################

# https://github.com/yfukasawa/LongQC
# rule ont_longQC_preTrim:
#     input:
#         MERGED_FQ = "{OUTDIR}/{SAMPLE}/ont/merged.fq.gz"
#     output:
#         DIR = directory("{OUTDIR}/{SAMPLE}/longqc/ont_preAdapterScan")
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
        R1_FQ="{OUTDIR}/{SAMPLE}/tmp/ont/adapter_scan_readids/merged_adapter_R1.fq.gz",
        # R1_FQ_Trimmed = "{OUTDIR}/{SAMPLE}/tmp/ont/adapter_scan_readids/full_len_internalTrim_R1.fq.gz",
        R2_FQ="{OUTDIR}/{SAMPLE}/tmp/ont/adapter_scan_readids/merged_adapter_R2.fq.gz",
    output:
        R1_FQ=temp("{OUTDIR}/{SAMPLE}/tmp/ont/cut_R1.fq.gz"),
        R2_FQ=temp("{OUTDIR}/{SAMPLE}/tmp/ont/cut_R2.fq.gz"),
        JSON="{OUTDIR}/{SAMPLE}/ont/cutadapt.json",
    params:
        ADAPTER=config["R1_INTERNAL_ADAPTER"],  # Curio R1 internal adapter
        R1=config["R1_PRIMER"],  # R1 PCR primer (Visium & Seeker)
        ADAPTER_COUNT=4,  # number of adapters that can be trimmed from each read
        QUALITY_MIN=5,  # Low Q score for ONT...
        MIN_R2_LENGTH=12,
        OVERLAP=5,
        HOMOPOLYMER_ERROR_RATE=0.2,  # default error rate is 0.1
        HETEROPOLYMER_ERROR_RATE=0.2,  # default error rate is 0.1
        POLYA="A" * 100,
        POLYT="T" * 100,
        TSO="AAGCTGGTATCAACGCAGAGTGAATGGG",  # SlideSeq TSO - remove any polyadenylated TSOs
        TXG_TSO="AAGCAGTGGTATCAACGCAGAGTACATGGG",  # 10x TSO - remove any polyadenylated TSOs
        rcTSO="CCCATTCACTCTGCGTTGATACCAGCTT",  # rev comp of SlideSeq TSO
        rcTXG_TSO="CCCATGTACTCTGCGTTGATACCACTGCTT",  # rev-comp of 10x TSO sequence
        SEEKER_BB_LINKER="TCTTCAGCGTTCCCGAGA",  # Adapter between BB1 & BB2 in R1 
        rcSEEKER_BB_ADAPTER="AGAGCCCTTGCGACTTCT",  # Reverse of the adapter between BB1 & BB2 in R1 
    # min([config['CORES'],8]) # 8 core max
    threads: config["CORES"]
    log:
        log="{OUTDIR}/{SAMPLE}/ont/cutadapt.log",
    run:
        R1_LENGTHS = [
            RECIPE_SHEET["R1.finalLength"][recipe]
            for recipe in RECIPE_ONT_DICT[wildcards.SAMPLE]
        ]
        # R1_LENGTH = RECIPE_SHEET["R1.finalLength"][recipe]
        R1_LENGTH = min(R1_LENGTHS)  # + len(params.R1) # Add R1 primer length for ONT

        R1 = input.R1_FQ
        # R1 = input.R1_FQ_Trimmed
        R2 = input.R2_FQ

        # TODO - make some params recipe-specific
        shell(
            f"""
            mkdir -p $(dirname  {output.R1_FQ})

            echo "Minimum R1 length: {R1_LENGTH}" > {log.log}
            echo "" >> {log.log}

            {EXEC['CUTADAPT']} \
                --minimum-length {R1_LENGTH}: \
                --error-rate {params.HETEROPOLYMER_ERROR_RATE} \
                --overlap {params.OVERLAP} \
                --match-read-wildcards \
                --times {params.ADAPTER_COUNT} \
                -g READ1=X{params.R1} \
                -A POLYA_3p="{params.POLYA}X;max_error_rate={params.HOMOPOLYMER_ERROR_RATE}" \
                -B TSO={params.TSO} \
                -B TXG_TSO={params.TXG_TSO} \
                -B SEEKER_LINKER={params.SEEKER_BB_LINKER} \
                --pair-filter=any \
                -o {output.R1_FQ} \
                -p {output.R2_FQ} \
                --cores {threads} \
                --json {output.JSON} \
                {R1} {R2} \
            1>> {log.log}
            """
        )
        # --minimum-length {R1_LENGTH}:{params.MIN_R2_LENGTH} \
        # --maximum-length {2*R1_LENGTH}: \
        # -A POLYT_3p="{params.POLYT};max_error_rate={params.HOMOPOLYMER_ERROR_RATE}" \

        # -G {params.FIVE_PRIME_R2_SEEKER_BB_ADAPTER} \
        # -B {params.rcTSO} \
        # -B {params.rcTXG_TSO} \
        # --quality-cutoff {params.QUALITY_MIN} \
        # -a "{params.POLYT};max_error_rate={params.HOMOPOLYMER_ERROR_RATE}" \ # 3' poly(T) trimming on R1?



# TODO- internal adapter trimming for ONT/R1/slide-seq


# Trimming for R1 to handle Curio adapter issues. See README for recipe details (#TODO)
rule ont_R1_trimming:
    input:
        R1_FQ="{OUTDIR}/{SAMPLE}/tmp/ont/cut_R1.fq.gz",
    output:
        R1_FQ="{OUTDIR}/{SAMPLE}/tmp/ont/cut_trim_R1.fq.gz",  #temp()
    params:
        CB1end=8,  #TODO- move to config? or recipe_sheet?
        CB2start=27,
        CB2end=42,
        INTERNAL_TRIM_QC_LOG="{OUTDIR}/{SAMPLE}/internal_trim_qc.txt",
        TMPDIR="{OUTDIR}/{SAMPLE}/tmp/seqkit",
        INTERNAL_ADAPTER=config["R1_INTERNAL_ADAPTER"],  # Curio R1 internal adapter
    threads: config["CORES"]
    log:
        log="{OUTDIR}/{SAMPLE}/ont/R1_trimming.log",
    run:
        recipe = RECIPE_DICT[wildcards.SAMPLE]
        R1_LENGTH = RECIPE_SHEET["R1.finalLength"][recipe]

        # param handling for different alignment strategies
        if "hardTrim" in recipe:
            shell(  # "Hard" trimming, to remove the adapter based on hard-coded base positions
                f"""
                zcat {input.R1_FQ} \
                | awk -v s={params.CB1end} \
                    -v S={params.CB2start} \
                    -v E={params.CB2end} \
                    -f scripts/awk/hardTrimFq.awk \
                > {output.R1_FQ.strip('.gz')}

                {EXEC['PIGZ']} -f -p{threads} {output.R1_FQ.strip('.gz')}

                echo "Hard trimming performed on {input.R1_FQ}" > {log.log}
                """
            )
        elif "internalTrim" in recipe:
            # TODO- rewrite/speed up internal trimming!
            shell(  # Internal trimming to cut out the SlideSeq adapter sequence
                f"""
                python scripts/py/internal_adapter_trim_R1.py \
                    {params.INTERNAL_ADAPTER} \
                    {params.INTERNAL_TRIM_QC_LOG} \
                    {threads} \
                    {params.TMPDIR} \
                    {R1} \
                    {output.R1_FQ} \
                | tee {log.log}
                """
            )
        else:
            shell(  # Rename R1_FQ if no trimming needed
                f"""
                cp {input.R1_FQ} {output.R1_FQ} 
                echo "No trimming performed on {input.R1_FQ}..." {log.log}
                """
            )
