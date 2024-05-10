#############################################
## Read trimming rules
#############################################


# TSO & homopolymer trimming
# TODO: add "{ADAPTER};noindels" to adapter sequence trimming? - *Note- do not do this for the BB_ADAPTER
rule cutadapt:
    input:
        R1_FQ="{OUTDIR}/{SAMPLE}/tmp/merged_R1.fq.gz",
        R2_FQ="{OUTDIR}/{SAMPLE}/tmp/merged_R2.fq.gz",
    output:
        R1_FQ=temp("{OUTDIR}/{SAMPLE}/tmp/cut_R1.fq.gz"),
        R2_FQ=temp("{OUTDIR}/{SAMPLE}/tmp/cut_R2.fq.gz"),
        JSON="{OUTDIR}/{SAMPLE}/misc_logs/cutadapt1.json",
    params:
        RECIPE=lambda w: get_recipes(w, mode="ILMN"),
        R1_LENGTHS=lambda w: get_recipe_info(w, info_col="R1.finalLength", mode="ILMN"),
        # R1_LENGTH = 50,
        QUALITY_MIN=20,
        MIN_R2_LENGTH=12,
        OVERLAP=5,
        HOMOPOLYMER_ERROR_RATE=0.2,  # default error rate is 0.1
        POLYA="A" * 100,
        POLYT="T" * 100,
        NEXTERA="CTGTCTCTTATA",  # Nextera sequence
        rcNEXTERA="TATAAGAGACAG",  # Rev Comp of Nextera sequence        
        TSO="AAGCTGGTATCAACGCAGAGTGAATGGG",  # SlideSeq TSO - remove any polyadenylated TSOs
        TXG_TSO="AAGCAGTGGTATCAACGCAGAGTACATGGG",  # 10x TSO - remove any polyadenylated TSOs
        rcTSO="CCCATTCACTCTGCGTTGATACCAGCTT",  # rev comp of SlideSeq TSO
        rcTXG_TSO="CCCATGTACTCTGCGTTGATACCACTGCTT",  # rev-comp of 10x TSO sequence        
        ILMN_UNIVERSAL="AGATCGGAAGAG",  # Illumina Universal
        FIVE_PRIME_R2_TSO="CCCATTCACTCTGCGTTGATACCAGCTT",  # rev comp of SlideSeq TSO
        FIVE_PRIME_R2_TXG_TSO="CCCATGTACTCTGCGTTGATACCACTGCTT",  # rev-comp of 10x TSO sequence
        SEEKER_BB_ADAPTER="TCTTCAGCGTTCCCGAGA",  # Adapter between BB1 & BB2 in R1 
        rcSEEKER_BB_ADAPTER="AGAGCCCTTGCGACTTCT",  # Reverse of the adapter between BB1 & BB2 in R1 
    threads: config["CORES"]
    log:
        log="{OUTDIR}/{SAMPLE}/misc_logs/cutadapt1.log",
    run:
        R1_LENGTH = min(params.R1_LENGTHS)

        shell(
            f"""
            {EXEC['CUTADAPT']} \
                --minimum-length {R1_LENGTH}:{params.MIN_R2_LENGTH} \
                --quality-cutoff {params.QUALITY_MIN} \
                --overlap {params.OVERLAP} \
                --match-read-wildcards \
                --nextseq-trim=20 \
                -A POLYA_3p="{params.POLYA}X;max_error_rate={params.HOMOPOLYMER_ERROR_RATE}" \
                -A POLYT_3p="{params.POLYT}X;max_error_rate={params.HOMOPOLYMER_ERROR_RATE}" \
                -B TSO={params.TSO} \
                -B TXG_TSO={params.TXG_TSO} \
                -B SEEKER_ADAPTER="{params.SEEKER_BB_ADAPTER}" \
                -B NEXTERA="{params.NEXTERA}" \
                -A ILMN_UNIVERSAL_3p="{params.ILMN_UNIVERSAL}" \
                --pair-filter=any \
                -o {output.R1_FQ} \
                -p {output.R2_FQ} \
                --cores {threads} \
                --json {output.JSON} \
                {input.R1_FQ} {input.R2_FQ} \
            1> {log.log}
            """
        )


rule cutadapt2:
    input:
        R1_FQ="{OUTDIR}/{SAMPLE}/tmp/cut_R1.fq.gz",
        R2_FQ="{OUTDIR}/{SAMPLE}/tmp/cut_R2.fq.gz",
    output:
        R1_FQ=temp("{OUTDIR}/{SAMPLE}/tmp/twiceCut_R1.fq.gz"),
        R2_FQ=temp("{OUTDIR}/{SAMPLE}/tmp/twiceCut_R2.fq.gz"),
        JSON="{OUTDIR}/{SAMPLE}/misc_logs/cutadapt2.json",
    params:
        RECIPE=lambda w: get_recipes(w, mode="ILMN"),
        R1_LENGTHS=lambda w: get_recipe_info(w, info_col="R1.finalLength", mode="ILMN"),
        # R1_LENGTH = 50,
        QUALITY_MIN=20,
        MIN_R2_LENGTH=12,
        OVERLAP=5,
        HOMOPOLYMER_ERROR_RATE=0.2,  # default error rate is 0.1
        POLYA="A" * 100,
        POLYT="T" * 100,
        NEXTERA="CTGTCTCTTATA",  # Nextera sequence
        rcNEXTERA="TATAAGAGACAG",  # Rev Comp of Nextera sequence        
        TSO="AAGCTGGTATCAACGCAGAGTGAATGGG",  # SlideSeq TSO - remove any polyadenylated TSOs
        TXG_TSO="AAGCAGTGGTATCAACGCAGAGTACATGGG",  # 10x TSO - remove any polyadenylated TSOs
        rcTSO="CCCATTCACTCTGCGTTGATACCAGCTT",  # rev comp of SlideSeq TSO
        rcTXG_TSO="CCCATGTACTCTGCGTTGATACCACTGCTT",  # rev-comp of 10x TSO sequence        
        ILMN_UNIVERSAL="AGATCGGAAGAG",  # Illumina Universal
        FIVE_PRIME_R2_TSO="CCCATTCACTCTGCGTTGATACCAGCTT",  # rev comp of SlideSeq TSO
        FIVE_PRIME_R2_TXG_TSO="CCCATGTACTCTGCGTTGATACCACTGCTT",  # rev-comp of 10x TSO sequence
        SEEKER_BB_ADAPTER="TCTTCAGCGTTCCCGAGA",  # Adapter between BB1 & BB2 in R1 
        rcSEEKER_BB_ADAPTER="AGAGCCCTTGCGACTTCT",  # Reverse of the adapter between BB1 & BB2 in R1 
    threads: config["CORES"]
    log:
        log="{OUTDIR}/{SAMPLE}/misc_logs/cutadapt2.log",
    run:
        R1_LENGTH = min(params.R1_LENGTHS)

        shell(
            f"""
            {EXEC['CUTADAPT']} \
                --minimum-length {R1_LENGTH}:{params.MIN_R2_LENGTH} \
                --quality-cutoff {params.QUALITY_MIN} \
                --overlap {params.OVERLAP} \
                --match-read-wildcards \
                --nextseq-trim=20 \
                -A POLYA_3p="{params.POLYA}X;max_error_rate={params.HOMOPOLYMER_ERROR_RATE}" \
                -A POLYT_3p="{params.POLYT}X;max_error_rate={params.HOMOPOLYMER_ERROR_RATE}" \
                -B TSO={params.TSO} \
                -B TXG_TSO={params.TXG_TSO} \
                -B SEEKER_ADAPTER="{params.SEEKER_BB_ADAPTER}" \
                -B NEXTERA="{params.NEXTERA}" \
                -A ILMN_UNIVERSAL_3p="{params.ILMN_UNIVERSAL}" \
                --pair-filter=any \
                -o {output.R1_FQ} \
                -p {output.R2_FQ} \
                --cores {threads} \
                --json {output.JSON} \
                {input.R1_FQ} {input.R2_FQ} \
            1> {log.log}
            """
        )


# Trimming for R1 to handle Curio adapter issues. See README for recipe details
## "Hard" trimming, to remove the adapter based on hard-coded base positions
rule R1_hardTrimming:
    input:
        R1_FQ="{OUTDIR}/{SAMPLE}/tmp/twiceCut_R1.fq.gz",
    output:
        R1_FQ="{OUTDIR}/{SAMPLE}/tmp/twiceCut_hardTrim_R1.fq.gz",  #temp()
    params:
        CB1end=8,  #TODO- move to config? or recipe_sheet?
        CB2start=27,
        CB2end=42,
        INTERNAL_TRIM_QC_LOG="{OUTDIR}/{SAMPLE}/internal_trim_qc.txt",
        TMPDIR="{OUTDIR}/{SAMPLE}/tmp/seqkit",
        RECIPE=lambda w: get_recipes(w, mode="ILMN"),
        R1_LENGTH=lambda w: get_recipe_info(w, info_col="R1.finalLength", mode="ILMN"),
    threads: config["CORES"]
    log:
        log="{OUTDIR}/{SAMPLE}/misc_logs/R1_hardTrimming.log",
    run:
        shell(
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


## Internal trimming to cut out adapter sequences
###TODO- rewrite/speed up internal trimming!
rule R1_internalTrimming:
    input:
        R1_FQ="{OUTDIR}/{SAMPLE}/tmp/twiceCut_R1.fq.gz",
    output:
        R1_FQ="{OUTDIR}/{SAMPLE}/tmp/twiceCut_internalTrim_R1.fq.gz",  #temp()
        INTERNAL_TRIM_QC_LOG="{OUTDIR}/{SAMPLE}/misc_logs/internal_trim_qc.txt",
    params:
        CB1end=8,  #TODO- move to config? or recipe_sheet?
        CB2start=27,
        CB2end=42,
        TMPDIR="{OUTDIR}/{SAMPLE}/tmp/seqkit",
        ADAPTER=lambda w: get_recipe_info(w, "internal.adapter")[0],
        RECIPE=lambda w: get_recipes(w, mode="ILMN"),
        R1_LENGTH=lambda w: get_recipe_info(w, info_col="R1.finalLength", mode="ILMN"),
    threads: config["CORES"]
    log:
        log="{OUTDIR}/{SAMPLE}/misc_logs/R1_internalTrimming.log",
    run:
        # TODO- rewrite/speed up internal trimming!
        shell(  # Internal trimming to cut out the SlideSeq adapter sequence
            f"""
            python scripts/py/fastq_internal_adapter_trim_R1.py \
                {params.ADAPTER} \
                {output.INTERNAL_TRIM_QC_LOG} \
                {threads} \
                {params.TMPDIR} \
                {input.R1_FQ} \
                {output.R1_FQ} \
            | tee {log.log}
            """
        )


# rule R1_trimming:
#     input:
#         R1_FQ="{OUTDIR}/{SAMPLE}/tmp/twiceCut_R1.fq.gz",
#     output:
#         R1_FQ="{OUTDIR}/{SAMPLE}/tmp/twiceCut_Trim_R1.fq.gz",  #temp()
#     params:
#         CB1end=8,  #TODO- move to config? or recipe_sheet?
#         CB2start=27,
#         CB2end=42,
#         INTERNAL_TRIM_QC_LOG="{OUTDIR}/{SAMPLE}/internal_trim_qc.txt",
#         TMPDIR="{OUTDIR}/{SAMPLE}/tmp/seqkit",
#         INTERNAL_ADAPTER=config["R1_INTERNAL_ADAPTER"],  # Curio R1 internal adapter
#         RECIPE = lambda w: get_recipes(w, mode="ILMN"),
#         R1_LENGTH = lambda w: get_recipe_info(w, info_col="R1.finalLength", mode="ILMN")
#     threads: config["CORES"]
#     log:
#         log="{OUTDIR}/{SAMPLE}/R1_trimming.log",
#     run:
#         # param handling for different alignment strategies
#         for recipe in params.RECIPE:
#             if "hardTrim" in recipe:
#                 shell(  # "Hard" trimming, to remove the adapter based on hard-coded base positions
#                     f"""
#                     zcat {input.R1_FQ} \
#                     | awk -v s={params.CB1end} \
#                         -v S={params.CB2start} \
#                         -v E={params.CB2end} \
#                         -f scripts/awk/hardTrimFq.awk \
#                     > {output.R1_FQ.strip('.gz')}
#                     {EXEC['PIGZ']} -f -p{threads} {output.R1_FQ.strip('.gz')}
#                     echo "Hard trimming performed on {input.R1_FQ}" > {log.log}
#                     """
#                 )
#             elif "internalTrim" in recipe:
#                 # TODO- rewrite/speed up internal trimming!
#                 shell(  # Internal trimming to cut out the SlideSeq adapter sequence
#                     f"""
#                     python scripts/py/internal_adapter_trim_R1.py \
#                         {params.INTERNAL_ADAPTER} \
#                         {params.INTERNAL_TRIM_QC_LOG} \
#                         {threads} \
#                         {params.TMPDIR} \
#                         {R1} \
#                         {output.R1_FQ} \
#                     | tee {log.log}
#                     """
#                 )
#             else:
#                 shell(  # Rename R1_FQ if no trimming needed
#                     # cp {input.R1_FQ} {output.R1_FQ}
#                     f"""
#                     touch {output.R1_FQ}
#                     echo "No trimming performed on {input.R1_FQ}..." {log.log}
#                     """
#                 )
