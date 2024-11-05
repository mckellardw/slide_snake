#############################################
## Read trimming rules
#############################################
# TODO - make some params recipe-specific
rule ont_1b_cutadapt:
    input:
        R1_FQ="{OUTDIR}/{SAMPLE}/ont/tmp/merged_adapter_R1.fq.gz",
        R2_FQ="{OUTDIR}/{SAMPLE}/ont/tmp/merged_adapter_R2.fq.gz",
    output:
        R1_FQ=temp("{OUTDIR}/{SAMPLE}/ont/tmp/cut_R1.fq.gz"),
        R2_FQ=temp("{OUTDIR}/{SAMPLE}/ont/tmp/cut_R2.fq.gz"),
        JSON="{OUTDIR}/{SAMPLE}/ont/misc_logs/1b_cutadapt.json",
    params:
        RECIPE=lambda w: get_recipes(w, mode="ONT"),
        BC_PRIMER=lambda w: get_recipe_info(w, "fwd_primer", mode="list")[0],
        R1_LENGTH=lambda w: min(
            get_recipe_info(w, info_col="R1_finalLength", mode="list")
        )
        + len(get_recipe_info(w, "fwd_primer", mode="list")[0]),
        ADAPTER=lambda w: get_recipe_info(w, "internal_adapter", mode="ONT"),
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
        uMRT_TSO="CCCTCTCTCTCTCTTTCCTCTCTC",  #ISOPCR sequence from uMRT protocol; does not contain the 4T overhang
        SEEKER_BB_LINKER="TCTTCAGCGTTCCCGAGA",  # Adapter between BB1 & BB2 in R1 
        rcSEEKER_BB_ADAPTER="AGAGCCCTTGCGACTTCT",  # Reverse of the adapter between BB1 & BB2 in R1
        VNP="ACTTGCCTGTCGCTCTATCTTCTTTTT",
        SSP="TTTCTGTTGGTGCTGATATTGCT",
    resources:
        mem="16G",
    threads: config["CORES"]
    log:
        log="{OUTDIR}/{SAMPLE}/ont/misc_logs/1b_cutadapt.log",
    conda:
        f"{workflow.basedir}/envs/cutadapt.yml"
    shell:
        """
        mkdir -p $(dirname  {output.R1_FQ})

        echo "Minimum R1 length: {params.R1_LENGTH}" > {log.log}
        echo "" >> {log.log}

        cutadapt \
            --minimum-length {params.R1_LENGTH}:{params.MIN_R2_LENGTH} \
            --error-rate {params.HETEROPOLYMER_ERROR_RATE} \
            --overlap {params.OVERLAP} \
            --match-read-wildcards \
            --times {params.ADAPTER_COUNT} \
            -A POLYA_3p="{params.POLYA}X;max_error_rate={params.HOMOPOLYMER_ERROR_RATE}" \
            -B TSO={params.TSO} \
            -B TXG_TSO={params.TXG_TSO} \
            -B uMRT_TSO={params.uMRT_TSO} \
            -B SEEKER_LINKER={params.SEEKER_BB_LINKER} \
            --pair-filter=any \
            -o {output.R1_FQ} \
            -p {output.R2_FQ} \
            --cores {threads} \
            --json {output.JSON} \
            {input.R1_FQ} {input.R2_FQ} \
        1>> {log.log}
        """
        # -g R1_PRIMER=X{params.R1_PRIMER} \
        # --minimum-length {R1_LENGTH}:{params.MIN_R2_LENGTH} \
        # --maximum-length {2*R1_LENGTH}: \
        # -A POLYT_3p="{params.POLYT};max_error_rate={params.HOMOPOLYMER_ERROR_RATE}" \


# Trimming for R1 to handle adapter issues. See README for recipe details (#TODO)
## "Hard" trimming, to remove the adapter based on hard-coded base positions
rule ont_1b_R1_hardTrimming:
    input:
        R1_FQ="{OUTDIR}/{SAMPLE}/ont/tmp/cut_R1.fq.gz",
    output:
        R1_FQ="{OUTDIR}/{SAMPLE}/ont/tmp/cut_hardTrim_R1.fq.gz",
    params:
        CB1end=8,  #TODO- move to config? or recipe_sheet?
        CB2start=27,
        CB2end=42,
    resources:
        mem="16G",
    threads: config["CORES"]
    log:
        log="{OUTDIR}/{SAMPLE}/ont/misc_logs/1b_R1_hardTrimming.log",
    shell:
        """
        zcat {input.R1_FQ} \
        | awk -v s={params.CB1end} \
            -v S={params.CB2start} \
            -v E={params.CB2end} \
            -f scripts/awk/hardTrimFq.awk \
        > {output.R1_FQ.strip('.gz')}

        pigz -f -p{threads} {output.R1_FQ.strip('.gz')}

        echo "Hard trimming performed on {input.R1_FQ}" > {log.log}
        """


## Internal trimming to cut out adapter sequences
rule ont_1b_R1_internalTrim:
    input:
        R1_FQ="{OUTDIR}/{SAMPLE}/ont/tmp/cut_R1.fq.gz",
    output:
        R1_FQ="{OUTDIR}/{SAMPLE}/ont/tmp/cut_internalTrim_R1.fq.gz",
        # INTERNAL_TRIM_QC_LOG="{OUTDIR}/{SAMPLE}/ont/misc_logs/1b_internal_trim_qc.txt",
    params:
        TMPDIR="{OUTDIR}/{SAMPLE}/tmp/seqkit",
        ADAPTER=lambda w: get_recipe_info(w, info_col="internal_adapter", mode="list")[
            0
        ],
        # RECIPE = lambda w: get_recipes(w, mode="ONT"),
        R1_LENGTH=lambda w: get_recipe_info(w, info_col="R1_finalLength"),
        BC1_LENGTH=lambda w: max_sum_of_entries(
            get_recipe_info(w, info_col="BC_length", mode="ILMN")
        ),
        MIN_ALIGN_SCORE=10,
    resources:
        mem="16G",
    threads: config["CORES"]
    log:
        log="{OUTDIR}/{SAMPLE}/ont/misc_logs/1b_R1_internalTrimming.log",
    conda:
        f"{workflow.basedir}/envs/parasail.yml"
    shell:
        """
        python scripts/py/fastq_internal_adapter_trim_R1_v2.py \
            --adapter_seq {params.ADAPTER} \
            --n_cores {threads} \
            --tmp_dir {params.TMPDIR} \
            --fq1_in {input.R1_FQ} \
            --fq1_out {output.R1_FQ} \
            --min_adapter_start_pos {params.BC1_LENGTH} \
            --min_align_score {params.MIN_ALIGN_SCORE} \
        | tee {log.log}
        """


# cutaddapt clean-up after internalTrimming on R1 (remove reads which didn't have adapter)
rule ont_1b_cutadapt_internalTrimming:
    input:
        R1_FQ="{OUTDIR}/{SAMPLE}/ont/tmp/cut_internalTrim_R1.fq.gz",
        R2_FQ="{OUTDIR}/{SAMPLE}/ont/tmp/cut_R2.fq.gz",
    output:
        R1_FQ="{OUTDIR}/{SAMPLE}/ont/tmp/twiceCut_internalTrim_R1.fq.gz",
        R2_FQ="{OUTDIR}/{SAMPLE}/ont/tmp/twiceCut_internalTrim_R2.fq.gz",
        JSON="{OUTDIR}/{SAMPLE}/ont/misc_logs/1b_cutadapt_internalTrim.json",
    params:
        RECIPE=lambda w: get_recipes(w, mode="ONT"),
        R1_LENGTH=lambda w: min(
            get_recipe_info(w, info_col="R1_finalLength", mode="list")
        )
        + len(params.R1),
        # Add R1 primer length for ONT,
        ADAPTER=lambda w: get_recipe_info(w, "internal_adapter", mode="ONT"),
    resources:
        mem="16G",
    threads: config["CORES"]
    log:
        log="{OUTDIR}/{SAMPLE}/ont/misc_logs/1b_cutadapt_internalTrim.log",
    conda:
        f"{workflow.basedir}/envs/cutadapt.yml"
    shell:
        """
        mkdir -p $(dirname  {output.R1_FQ})

        echo "Minimum R1 length: {params.R1_LENGTH}" > {log.log}
        echo "" >> {log.log}

        cutadapt \
            --minimum-length {params.R1_LENGTH}: \
            --pair-filter=any \
            -o {output.R1_FQ} \
            -p {output.R2_FQ} \
            --cores {threads} \
            --json {output.JSON} \
            {input.R1_FQ} {input.R2_FQ} \
        1>> {log.log}
        """
