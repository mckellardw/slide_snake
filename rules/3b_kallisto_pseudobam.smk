# TODO- update for different recipes
rule kallisto_quant_bulk:
    input:
        R1_FQ_FILTERED = '{OUTDIR}/{sample}/tmp/final_filtered_R1.fq.gz',
        R2_FQ_FILTERED = '{OUTDIR}/{sample}/tmp/final_filtered_R2.fq.gz',
    output:
        GENOMEBAM = '{OUTDIR}/{sample}/kb/{RECIPE}/quant/pseudoalignments.bam'
    params:
        MEMLIMIT = config['MEMLIMIT']
    log:
        '{OUTDIR}/{sample}/kb/{RECIPE}/kallisto_quant_standard.log'
    threads:
        config['CORES']
    run:
        # tmp_recipe = RECIPE_DICT[wildcards.sample]
        KB_IDX = IDX_DICT[wildcards.sample]

        CHROMOSOMES = f"{REF_DICT[wildcards.sample]}/chrNameLength.txt"
        GENES_GTF = GTF_DICT[wildcards.sample]

        shell(
            f"""
            mkdir -p {OUTDIR}/{wildcards.sample}/kb/{RECIPE}/quant

            {EXEC['KALLISTO']} quant \
                --index {KB_IDX} \
                --output-dir {OUTDIR}/{wildcards.sample}/kb/{RECIPE}/quant/ \
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
