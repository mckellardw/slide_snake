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
# TODO- refactor to take info from recipe_sheet on barcode positions/lengths
# TODO- refactor to take a variable number of barcodes
# TODO- export to python script
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
        BC_LENGTHS=lambda w: get_recipe_info(w, info_col="BC_length", mode="list"),
    # resources:
    threads: 1
    run:
        recipes = "".join(get_recipes(wildcards, mode="list"))

        # load bc
        bc_map = pd.read_csv(input.BC_MAP, sep="\t", header=None)

        # split into multiple whitelists for separated barcodes
        # TODO- fix this code... Needs to be abstracted!
        if "seeker" in recipes:
            bc_1 = [bc[:8] for bc in list(bc_map[0])]
            bc_2 = [bc[8:] for bc in list(bc_map[0])]
            bc_us = [bc[:8] + "_" + bc[8:] for bc in list(bc_map[0])]

            bc_us_map = bc_map.copy()
            bc_us_map.iloc[:, 0] = bc_us

            # save bc files in {SAMPLE}/bc
            bc_us_map.to_csv(output.BC_US_MAP, header=False, index=False, sep="\t")
            pd.Series(bc_1).to_csv(output.BC_1, header=False, index=False)
            pd.Series(bc_2).to_csv(output.BC_2, header=False, index=False)
            pd.Series(list(set(bc_1))).to_csv(
                output.BC_UNIQ_1, header=False, index=False
            )
            pd.Series(list(set(bc_2))).to_csv(
                output.BC_UNIQ_2, header=False, index=False
            )
            pd.Series(bc_us).to_csv(output.BC_US, header=False, index=False)
        elif "decoder" in recipes:
            bc_1 = [bc[:8] for bc in list(bc_map[0])]
            bc_2 = [bc[8:] for bc in list(bc_map[0])]
            bc_us = [bc[:8] + "_" + bc[8:] for bc in list(bc_map[0])]

            bc_us_map = bc_map.copy()
            bc_us_map.iloc[:, 0] = bc_us

            # save bc files in {SAMPLE}/bc
            bc_us_map.to_csv(output.BC_US_MAP, header=False, index=False, sep="\t")
            pd.Series(bc_1).to_csv(output.BC_1, header=False, index=False)
            pd.Series(bc_2).to_csv(output.BC_2, header=False, index=False)
            pd.Series(list(set(bc_1))).to_csv(
                output.BC_UNIQ_1, header=False, index=False
            )
            pd.Series(list(set(bc_2))).to_csv(
                output.BC_UNIQ_2, header=False, index=False
            )
            pd.Series(bc_us).to_csv(output.BC_US, header=False, index=False)

        elif "microST" in recipes:
            bc_1 = [bc[:10] for bc in list(bc_map[0])]
            bc_2 = [bc[10:] for bc in list(bc_map[0])]
            bc_us = [bc[:10] + "_" + bc[10:] for bc in list(bc_map[0])]

            bc_us_map = bc_map.copy()
            bc_us_map.iloc[:, 0] = bc_us

            # save bc files in {SAMPLE}/bc
            bc_us_map.to_csv(output.BC_US_MAP, header=False, index=False, sep="\t")
            pd.Series(bc_1).to_csv(output.BC_1, header=False, index=False)
            pd.Series(bc_2).to_csv(output.BC_2, header=False, index=False)
            pd.Series(list(set(bc_1))).to_csv(
                output.BC_UNIQ_1, header=False, index=False
            )
            pd.Series(list(set(bc_2))).to_csv(
                output.BC_UNIQ_2, header=False, index=False
            )
            pd.Series(bc_us).to_csv(output.BC_US, header=False, index=False)

        else:
            shell(
                f"""
                cat {input.BC_MAP} > {output.BC_US_MAP}
                touch {output.BC_1}
                touch {output.BC_2}
                touch {output.BC_UNIQ_1}
                touch {output.BC_UNIQ_2}
                cat {input.BC_MAP} | cut -f1 > {output.BC_US}
                """
            )


# For Seeker - Insert the adapter sequence into the bead barcodes for barcode matching/alignment
rule BC_insert_adapter_into_list:
    input:
        BC_MAP="{OUTDIR}/{SAMPLE}/bc/map.txt",
    output:
        BC_ADAPTER="{OUTDIR}/{SAMPLE}/bc/whitelist_adapter.txt",
        BC_ADAPTER_R1="{OUTDIR}/{SAMPLE}/bc/whitelist_adapter_r1.txt",
        BC_ADAPTER_MAP="{OUTDIR}/{SAMPLE}/bc/map_adapter.txt",
        BC_ADAPTER_R1_MAP="{OUTDIR}/{SAMPLE}/bc/map_adapter_R1.txt",
    params:
        ADAPTER=lambda w: get_recipe_info(w, "internal_adapter", mode="list")[0],
        BC_PRIMER=lambda w: get_recipe_info(w, "fwd_primer", mode="list")[0],
        recipes_to_split=["seeker", "microST", "decoder"],
    # resources:
    threads: 1
    run:
        recipes = "".join(get_recipes(wildcards, mode="list"))

        if any(recipe in recipes for recipe in params.recipes_to_split):
            # load bc
            bc_map = pd.read_csv(input.BC_MAP, sep="\t", header=None)
            if "seeker" in recipes:
                # split for 2 separate barcodes
                bc_1 = [bc[:8] for bc in list(bc_map[0])]
                bc_2 = [
                    bc[8:] for bc in list(bc_map[0])
                ]  # Stitch bc_1, adapter, and bc_2

                bc_adapter = [
                    f"{item1}{params.ADAPTER}{item2}"
                    for item1, item2 in zip(bc_1, bc_2)
                ]
                bc_adapter_r1 = [
                    f"{params.BC_PRIMER}{item1}{params.ADAPTER}{item2}"
                    for item1, item2 in zip(bc_1, bc_2)
                ]

            elif "microST" in recipes:
                # split for 2 separate barcodes
                bc_1 = [bc[0:10] for bc in list(bc_map[0])]
                bc_2 = [bc[10:] for bc in list(bc_map[0])]

                bc_adapter = [
                    f"{item1}{params.ADAPTER}{item2}"
                    for item1, item2 in zip(bc_1, bc_2)
                ]
                bc_adapter_r1 = [
                    f"{params.BC_PRIMER}{item1}{params.ADAPTER}{item2}"
                    for item1, item2 in zip(bc_1, bc_2)
                ]

            elif "decoder" in recipes:
                print("TODO")

            bc_adapter_map = pd.DataFrame(list(zip(bc_adapter, bc_map[1], bc_map[2])))
            bc_adapter_r1_map = pd.DataFrame(
                list(zip(bc_adapter_r1, bc_map[1], bc_map[2]))
            )

            # save bc files in {SAMPLE}/bc
            pd.Series(bc_adapter).to_csv(
                output.BC_ADAPTER, sep="\n", header=False, index=False
            )
            pd.Series(bc_adapter_r1).to_csv(
                output.BC_ADAPTER_R1, sep="\n", header=False, index=False
            )

            bc_adapter_map.to_csv(
                output.BC_ADAPTER_MAP, sep="\t", header=False, index=False
            )
            bc_adapter_r1_map.to_csv(
                output.BC_ADAPTER_R1_MAP, sep="\t", header=False, index=False
            )
        else:
            shell(f"touch {' '.join(output)}")
