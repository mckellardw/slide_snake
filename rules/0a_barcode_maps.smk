# Rules for formatting barcode whitelists/maps for different chemistries


# Copy barcode map
rule BC_copy_barcode_map:
    input:
        BC_MAP=lambda wildcards: SAMPLE_SHEET["BC_map"][wildcards.SAMPLE],
    output:
        BC_MAP="{OUTDIR}/{SAMPLE}/bc/map.txt",
    # resources:
    threads: 1
    run:
        if input.BC_MAP.endswith(".gz"):
            shell(
                f"""
                mkdir -p $(dirname {output.BC_MAP})
                zcat {input.BC_MAP} > {output.BC_MAP}
                """
            )
        else:
            shell(
                f"""
                mkdir -p $(dirname {output.BC_MAP})
                cat {input.BC_MAP} > {output.BC_MAP}
                """
            )


rule BC_get_simple_whitelist:
    input:
        BC_MAP="{OUTDIR}/{SAMPLE}/bc/map.txt",
    output:
        BC="{OUTDIR}/{SAMPLE}/bc/whitelist.txt",
    # resources:
    threads: 1
    shell:
        """
        mkdir -p $(dirname {output.BC})
        cut -f1 {input.BC_MAP} > {output.BC}
        """


# Split the barcodes and save whitelists
rule BC_write_whitelist_variants:
    input:
        BC_MAP="{OUTDIR}/{SAMPLE}/bc/map.txt",
    output:
        BC_US_MAP="{OUTDIR}/{SAMPLE}/bc/map_underscore.txt",  # Barcode Map With Underscore for STAR
        BC_1="{OUTDIR}/{SAMPLE}/bc/whitelist_1.txt",  # Barcode #1
        BC_2="{OUTDIR}/{SAMPLE}/bc/whitelist_2.txt",  # Barcode #2
        BC_UNIQ_1="{OUTDIR}/{SAMPLE}/bc/whitelist_uniq_1.txt",  # barcode #1, unique values
        BC_UNIQ_2="{OUTDIR}/{SAMPLE}/bc/whitelist_uniq_2.txt",  # Barcode #2, unique values
        BC_US="{OUTDIR}/{SAMPLE}/bc/whitelist_underscore.txt",  # Barcode With Underscore for STAR
    params:
        BC_LENGTHS=lambda w: get_recipe_info(w, info_col="BC_length", mode="list")[0],
        RECIPES=lambda w: get_recipes(w, mode="list"),
    threads: 1
    log:
        log="{OUTDIR}/{SAMPLE}/bc/process_barcode_whitelist.log",
        err="{OUTDIR}/{SAMPLE}/bc/process_barcode_whitelist.err",
    shell:
        """
        python scripts/py/process_barcode_whitelist.py \
            --bc-map-file {input.BC_MAP} \
            --bc-us-map {output.BC_US_MAP} \
            --bc-1 {output.BC_1} \
            --bc-2 {output.BC_2} \
            --bc-uniq-1 {output.BC_UNIQ_1} \
            --bc-uniq-2 {output.BC_UNIQ_2} \
            --bc-us {output.BC_US} \
            --bc-lengths "{params.BC_LENGTHS}" \
        > {log.log} \
        2> {log.err}
        """
