# TODO- update for different recipes
rule kallisto_quant_bulk:
    input:
        FINAL_R1_FQ = '{OUTDIR}/{sample}/tmp/{sample}_R1_final_trimmed.fq.gz',
        FINAL_R2_FQ = '{OUTDIR}/{sample}/tmp/{sample}_R2_final_trimmed.fq.gz',
        BB = "{OUTDIR}/{sample}/tmp/whitelist.txt"
    output:
        GENOMEBAM = '{OUTDIR}/{sample}/kb/quant/pseudoalignments.bam'
    params:
        OUTDIR = config['OUTDIR'],
        KALLISTO_EXEC = config['KALLISTO_EXEC'],
        # KB_IDX = config['KB_IDX'],
        MEMLIMIT = config['MEMLIMIT']
    log:
        '{OUTDIR}/{sample}/kb/kallisto_quant_standard.log'
    threads:
        config['CORES']
    run:
        tmp_recipe = RECIPE_DICT[wildcards.sample]
        KB_IDX = IDX_DICT[wildcards.sample]
        BB_WHITELIST = f"{input.BB}" 

        STAR_REF = REF_DICT[wildcards.sample] + "/chrNameLength.txt"
        GENES_GTF = GTF_DICT[wildcards.sample]

        shell(
            f"""
            mkdir -p {params.OUTDIR}/{wildcards.sample}/kb/quant

            {params.KALLISTO_EXEC} quant \
            -i {KB_IDX} \
            -o {params.OUTDIR}/{wildcards.sample}/kb/quant/ \
            -t {threads} \
            --fr-stranded \
            --single \
            -l 85 \
            -s 10 \
            --genomebam \
            --chromosomes {params.CHROMOSOMES} \
            --gtf {GENES_GTF} \
            {input.FINAL_R2_FQ}
            """
        )
