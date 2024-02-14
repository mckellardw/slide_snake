# TODO- update for different recipes
rule kallisto_quant_bulk:
    input:
        R1_FQ = '{OUTDIR}/{SAMPLE}/tmp/cut_R1.fq.gz',
        R2_FQ = '{OUTDIR}/{SAMPLE}/tmp/cut_R2.fq.gz',
        R1_FQ_TWICE_CUT = '{OUTDIR}/{SAMPLE}/tmp/twiceCut_R1.fq.gz',
        R2_FQ_TWICE_CUT = '{OUTDIR}/{SAMPLE}/tmp/twiceCut_R2.fq.gz',
        R1_FQ_STAR_FILTERED = '{OUTDIR}/{SAMPLE}/rRNA/STARsolo/final_filtered_R1.fq.gz',
        R2_FQ_STAR_FILTERED = '{OUTDIR}/{SAMPLE}/rRNA/STARsolo/final_filtered_R2.fq.gz',
        R1_FQ_BWA_FILTERED  = '{OUTDIR}/{SAMPLE}/rRNA/bwa/final_filtered_R1.fq.gz',
        R2_FQ_BWA_FILTERED  = '{OUTDIR}/{SAMPLE}/rRNA/bwa/final_filtered_R2.fq.gz',
    output:
        GENOMEBAM = '{OUTDIR}/{SAMPLE}/kb/{RECIPE}/quant/pseudoalignments.bam'
    params:
        MEMLIMIT = config['MEMLIMIT']
    log:
        log = '{OUTDIR}/{SAMPLE}/kb/{RECIPE}/kallisto_quant.log'
    threads:
        config['CORES']
    run:
        # recipe = RECIPE_DICT[wildcards.SAMPLE]
        recipe = wildcards.RECIPE
        
        KB_IDX = IDX_DICT[wildcards.SAMPLE]        
        KB_X = RECIPE_SHEET["kb.x"][recipe]
        
        # Select input reads based on alignment recipe
        if "rRNA.STAR" in recipe: # Use trimmed & STAR-rRNA-filtered .fq's
            # R1 = input.R1_FQ_STAR_FILTERED
            R2 = input.R2_FQ_STAR_FILTERED
        elif "rRNA.bwa" in recipe: #TODO Use trimmed & bwa-rRNA-filtered .fq's
            # R1 = input.R1_FQ_BWA_FILTERED
            R2 = input.R2_FQ_BWA_FILTERED
        elif "rRNA" not in recipe: # just trimmed .fq's
            # R1 = input.R1_FQ
            # R2 = input.R2_FQ
            # R1 = input.R1_FQ_TWICE_CUT
            R2 = input.R2_FQ_TWICE_CUT
        else:
            print("I just don't know what to do with myself...")

        CHROMOSOMES = f"{REF_DICT[wildcards.SAMPLE]}/chrNameLength.txt"
        GENES_GTF = GTF_DICT[wildcards.SAMPLE]

        shell(
            f"""
            mkdir -p $(dirname {output.GENOMEBAM})

            {EXEC['KALLISTO']} quant \
                --index {KB_IDX} \
                --output-dir $(dirname {output.GENOMEBAM}) \
                --threads {threads} \
                --fr-stranded \
                --single \
                --fragment-length 85 \
                --sd 10 \
                --genomebam \
                --chromosomes {CHROMOSOMES} \
                --gtf {GENES_GTF} \
                {R2} \
            1> {log.log}
            """
        )
