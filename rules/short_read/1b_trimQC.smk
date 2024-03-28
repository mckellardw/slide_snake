#############################################
## Read trimming & QC rules
#############################################
# fastqc on R1 & R2 before trimming
rule fastQC_preTrim:
    input:
        MERGED_FQ="{OUTDIR}/{SAMPLE}/tmp/merged_{READ}.fq.gz",
    output:
        fastqcDir=directory("{OUTDIR}/{SAMPLE}/fastqc/preCutadapt_{READ}"),
    params:
        adapters=config["FASTQC_ADAPTERS"],
    threads: config["CORES"]
        # min([config['CORES'],8]) # 8 core max
    run:
        shell(
            f"""
            mkdir -p {output.fastqcDir}

            {EXEC['FASTQC']} \
                --outdir {output.fastqcDir} \
                --threads {threads} \
                -a {params.adapters} \
                {input.MERGED_FQ}
            """
        )


# Trimming for R1 to handle Curio adapter issues. See README for recipe details (#TODO)
rule R1_trimming:
    input:
        R1_FQ="{OUTDIR}/{SAMPLE}/tmp/merged_R1.fq.gz",
    output:
        R1_FQ="{OUTDIR}/{SAMPLE}/tmp/merged_trimmed_R1.fq.gz",  #temp()
    params:
        CB1end=8,  #TODO- move to config? or recipe_sheet?
        CB2start=27,
        CB2end=42,
        INTERNAL_TRIM_QC_LOG="{OUTDIR}/{SAMPLE}/internal_trim_qc.txt",
        TMPDIR="{OUTDIR}/{SAMPLE}/tmp/seqkit",
        INTERNAL_ADAPTER=config["R1_INTERNAL_ADAPTER"],  # Curio R1 internal adapter
    threads: config["CORES"]
    log:
        log="{OUTDIR}/{SAMPLE}/R1_trimming.log",
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


#


# TSO & homopolymer trimming
# TODO: add "{ADAPTER};noindels" to adapter sequence trimming? - *Note- do not do this for the BB_ADAPTER
rule cutadapt:
    input:
        R1_FQ="{OUTDIR}/{SAMPLE}/tmp/merged_R1.fq.gz",
        # R1_FQ_HardTrim = '{OUTDIR}/{SAMPLE}/tmp/merged_R1_HardTrim.fq.gz',
        # R1_FQ_InternalTrim = '{OUTDIR}/{SAMPLE}/tmp/merged_R1_InternalTrim.fq.gz',
        R1_FQ_Trimmed="{OUTDIR}/{SAMPLE}/tmp/merged_trimmed_R1.fq.gz",
        R2_FQ="{OUTDIR}/{SAMPLE}/tmp/merged_R2.fq.gz",
    output:
        FINAL_R1_FQ=temp("{OUTDIR}/{SAMPLE}/tmp/cut_R1.fq.gz"),
        FINAL_R2_FQ=temp("{OUTDIR}/{SAMPLE}/tmp/cut_R2.fq.gz"),
        JSON="{OUTDIR}/{SAMPLE}/cutadapt.json",
    params:
        # R1_LENGTH = 50,
        ADAPTER_COUNT=4,  # number of adapters that can be trimmed from each read
        QUALITY_MIN=20,
        MIN_R2_LENGTH=12,
        OVERLAP=5,
        HOMOPOLYMER_ERROR_RATE=0.2,  # default error rate is 0.1
        # THREE_PRIME_R1_POLYA = "A"*100,
        THREE_PRIME_R2_POLYA="A" * 100,
        # THREE_PRIME_R2_POLYG = "G"*100,
        THREE_PRIME_R2_POLYT="T" * 100,
        THREE_PRIME_R2_NEXTERA="CTGTCTCTTATA",  # Nextera sequence
        THREE_PRIME_R2_rcNEXTERA="TATAAGAGACAG",  # Rev Comp of Nextera sequence
        THREE_PRIME_R2_TSO="AAGCTGGTATCAACGCAGAGTGAATGGG",  # SlideSeq TSO - remove any polyadenylated TSOs
        THREE_PRIME_R2_TXG_TSO="AAGCAGTGGTATCAACGCAGAGTACATGGG",  # 10x TSO - remove any polyadenylated TSOs
        THREE_PRIME_R2_ILLUMINA_UNI="AGATCGGAAGAG",  # Illumina Universal
        FIVE_PRIME_R2_TSO="CCCATTCACTCTGCGTTGATACCAGCTT",  # rev comp of SlideSeq TSO
        FIVE_PRIME_R2_TXG_TSO="CCCATGTACTCTGCGTTGATACCACTGCTT",  # rev-comp of 10x TSO sequence
        THREE_PRIME_R2_SEEKER_BB_ADAPTER="TCTTCAGCGTTCCCGAGA",  # Adapter between BB1 & BB2 in R1 
        FIVE_PRIME_R2_SEEKER_BB_ADAPTER="AGAGCCCTTGCGACTTCT",  # Reverse of the adapter between BB1 & BB2 in R1 
    # min([config['CORES'],8]) # 8 core max
    threads: config["CORES"]
    log:
        log="{OUTDIR}/{SAMPLE}/cutadapt.log",
    run:
        recipe = str(
            RECIPE_DICT[wildcards.SAMPLE][0]
        )  # TODO- fix multi-R1 trimming handling

        R1_LENGTH = RECIPE_SHEET["R1.finalLength"][recipe]

        # TODO- fix multi-R1 trimming handling
        # param handling for different alignment strategies
        # if "Trim" in recipe:
        #     R1 = input.R1_FQ_Trimmed
        # else:
        #     R1 = input.R1_FQ

        R1 = input.R1_FQ_Trimmed
        R2 = input.R2_FQ

        shell(
            f"""
            {EXEC['CUTADAPT']} \
                --minimum-length {R1_LENGTH}:{params.MIN_R2_LENGTH} \
                --quality-cutoff {params.QUALITY_MIN} \
                --overlap {params.OVERLAP} \
                --match-read-wildcards \
                --nextseq-trim=20 \
                --times {params.ADAPTER_COUNT} \
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
                --json {output.JSON} \
                {R1} {R2} \
            1> {log.log}
            """
        )


#


# fastqc on R1 after linker removal & R2 trimming/filtering
rule fastQC_postTrim:
    input:
        FINAL_FQ="{OUTDIR}/{SAMPLE}/tmp/cut_{READ}.fq.gz",
    output:
        fastqcDir=directory("{OUTDIR}/{SAMPLE}/fastqc/postCutadapt_{READ}"),
        # fastqcReport = ''
    threads: config["CORES"]
        # min([config['CORES'],8]) # 8 core max
    params:
        adapters=config["FASTQC_ADAPTERS"],
    run:
        shell(
            f"""
            mkdir -p {output.fastqcDir}

            {EXEC['FASTQC']} \
                --outdir {output.fastqcDir} \
                --threads {threads} \
                -a {params.adapters} \
                {input.FINAL_FQ}
            """
        )


##########################################


rule cutadapt2:
    input:
        FINAL_R1_FQ="{OUTDIR}/{SAMPLE}/tmp/cut_R1.fq.gz",
        FINAL_R2_FQ="{OUTDIR}/{SAMPLE}/tmp/cut_R2.fq.gz",
    output:
        FINAL_R1_FQ=temp("{OUTDIR}/{SAMPLE}/tmp/twiceCut_R1.fq.gz"),
        FINAL_R2_FQ=temp("{OUTDIR}/{SAMPLE}/tmp/twiceCut_R2.fq.gz"),
        JSON="{OUTDIR}/{SAMPLE}/cutadapt2.json",
    params:
        # R1_LENGTH = 50,
        QUALITY_MIN=20,
        MIN_R2_LENGTH=12,
        OVERLAP=5,
        HOMOPOLYMER_ERROR_RATE=0.2,  # default error rate is 0.1
        # THREE_PRIME_R1_POLYA = "A"*100,
        THREE_PRIME_R2_POLYA="A" * 100,
        # THREE_PRIME_R2_POLYG = "G"*100,
        THREE_PRIME_R2_POLYT="T" * 100,
        THREE_PRIME_R2_NEXTERA="CTGTCTCTTATA",  # Nextera sequence
        THREE_PRIME_R2_rcNEXTERA="TATAAGAGACAG",  # Rev Comp of Nextera sequence
        THREE_PRIME_R2_TSO="AAGCTGGTATCAACGCAGAGTGAATGGG",  # SlideSeq TSO - remove any polyadenylated TSOs
        THREE_PRIME_R2_TXG_TSO="AAGCAGTGGTATCAACGCAGAGTACATGGG",  # 10x TSO - remove any polyadenylated TSOs
        THREE_PRIME_R2_ILLUMINA_UNI="AGATCGGAAGAG",  # Illumina Universal
        FIVE_PRIME_R2_TSO="CCCATTCACTCTGCGTTGATACCAGCTT",  # rev comp of SlideSeq TSO
        FIVE_PRIME_R2_TXG_TSO="CCCATGTACTCTGCGTTGATACCACTGCTT",  # rev-comp of 10x TSO sequence
        THREE_PRIME_R2_SEEKER_BB_ADAPTER="TCTTCAGCGTTCCCGAGA",  # Adapter between BB1 & BB2 in R1 
        FIVE_PRIME_R2_SEEKER_BB_ADAPTER="AGAGCCCTTGCGACTTCT",  # Reverse of the adapter between BB1 & BB2 in R1 
    # min([config['CORES'],8]) # 8 core max
    threads: config["CORES"]
    log:
        log="{OUTDIR}/{SAMPLE}/cutadapt_round2.log",
    run:
        recipe = str(
            RECIPE_DICT[wildcards.SAMPLE][0]
        )  # TODO- fix multi-R1 trimming handling

        R1_LENGTH = RECIPE_SHEET["R1.finalLength"][recipe]

        R1 = input.FINAL_R1_FQ
        R2 = input.FINAL_R2_FQ

        shell(
            f"""
            {EXEC['CUTADAPT']} \
                --minimum-length {R1_LENGTH}:{params.MIN_R2_LENGTH} \
                --quality-cutoff {params.QUALITY_MIN} \
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
                --json {output.JSON} \
                {R1} {R2} \
            1> {log.log}
            """
        )


rule fastQC_twiceTrim:
    input:
        FINAL_FQ="{OUTDIR}/{SAMPLE}/tmp/twiceCut_{READ}.fq.gz",
    output:
        fastqcDir=directory("{OUTDIR}/{SAMPLE}/fastqc/twiceCutadapt_{READ}"),
        # fastqcReport = ''
    threads: config["CORES"]
        # min([config['CORES'],8]) # 8 core max
    params:
        adapters=config["FASTQC_ADAPTERS"],
    run:
        shell(
            f"""
            mkdir -p {output.fastqcDir}

            {EXEC['FASTQC']} \
                --outdir {output.fastqcDir} \
                --threads {threads} \
                -a {params.adapters} \
                {input.FINAL_FQ}
            """
        )
