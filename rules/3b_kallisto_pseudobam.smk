# TODO- update for different recipes
rule kallisto_quant_bulk:
    input:
        FINAL_R1_FQ = '{OUTDIR}/{sample}/tmp/{sample}_R1_final_trimmed.fq.gz',
        FINAL_R2_FQ = '{OUTDIR}/{sample}/tmp/{sample}_R2_final_trimmed.fq.gz'
    output:
        GENOMEBAM = '{OUTDIR}/{sample}/kb/quant/pseudoalignments.bam'
    params:
        MEMLIMIT = config['MEMLIMIT']
    log:
        '{OUTDIR}/{sample}/kb/kallisto_quant_standard.log'
    threads:
        config['CORES']
    run:
        # tmp_recipe = RECIPE_DICT[wildcards.sample]
        KB_IDX = IDX_DICT[wildcards.sample]

        CHROMOSOMES = f"{REF_DICT[wildcards.sample]}/chrNameLength.txt"
        GENES_GTF = GTF_DICT[wildcards.sample]

        shell(
            f"""
            mkdir -p {OUTDIR}/{wildcards.sample}/kb/quant

            {KALLISTO_EXEC} quant \
            --index {KB_IDX} \
            --output-dir {OUTDIR}/{wildcards.sample}/kb/quant/ \
            --threads {threads} \
            --fr-stranded \
            --single \
            --fragment-length 85 \
            --sd 10 \
            --genomebam \
            --chromosomes {CHROMOSOMES} \
            --gtf {GENES_GTF} \
            {input.FINAL_R2_FQ}
            """
        )
