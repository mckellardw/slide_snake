# TODO- update for different recipes
rule kallisto_quant_bulk:
    input:
        R1_FQ_FILTERED = '{OUTDIR}/{SAMPLE}/tmp/final_filtered_R1.fq.gz',
        R2_FQ_FILTERED = '{OUTDIR}/{SAMPLE}/tmp/final_filtered_R2.fq.gz',
    output:
        GENOMEBAM = '{OUTDIR}/{SAMPLE}/kb/{RECIPE}/quant/pseudoalignments.bam'
    params:
        MEMLIMIT = config['MEMLIMIT']
    log:
        log = '{OUTDIR}/{SAMPLE}/kb/{RECIPE}/kallisto_quant.log'
    threads:
        config['CORES']
    run:
        # tmp_recipe = RECIPE_DICT[wildcards.SAMPLE]
        KB_IDX = IDX_DICT[wildcards.SAMPLE]

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
                {input.FINAL_R2_FQ} \
            1> {log.log}
            """
        )
