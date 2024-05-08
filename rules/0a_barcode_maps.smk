# Rules for formatting barcode whitelists/maps for different chemistries


# Copy barcode map
rule copy_barcode_map:
    input:
        BC_MAP=lambda wildcards: BC_DICT[wildcards.SAMPLE],
    output:
        BC_MAP="{OUTDIR}/{SAMPLE}/bb/whitelist_map.txt",
    threads: 1
    run:
        shell(f"cp {input.BC_MAP} {output.BC_MAP}")


rule get_simple_whitelist:
    input:
        BC_MAP="{OUTDIR}/{SAMPLE}/bb/whitelist_map.txt",
    output:
        BC="{OUTDIR}/{SAMPLE}/bb/whitelist.txt",
    threads: 1
    run:
        shell(f"cut -f1 {input.BC_MAP} > {output.BC}")


# Split the barcodes and save whitelists
# TODO add function for handling variable whitelist outputs...
# TODO- refactor to take info from recipe_sheet on barcode positions/lengths
rule write_whitelist_variants:
    input:
        BC_MAP="{OUTDIR}/{SAMPLE}/bb/whitelist_map.txt",
    output:
        BC_1="{OUTDIR}/{SAMPLE}/bb/whitelist_1.txt",  # Bead barcode #1
        BC_2="{OUTDIR}/{SAMPLE}/bb/whitelist_2.txt",  # Bead Barcode #2
        BC_US="{OUTDIR}/{SAMPLE}/bb/whitelist_underscore.txt",  # Bead Barcode With Underscore for STAR
    threads: 1
    run:
        recipes = get_recipes(wildcards, mode="concatenate")[0]
        print(recipes)

        # load bb
        bc_df = pd.read_csv(input.BC_MAP, sep="\t", header=None).iloc[:, 0]

        # split into multiple whitelists for separated barcodes
        if "seeker" in recipes or "decoder" in recipes:
            bc_1 = [bb[:8] for bb in list(bc_df.values)]
            bc_1 = list(set(bc_1))

            bc_2 = [bb[8:] for bb in list(bc_df.values)]
            bc_2 = list(set(bc_2))

            bc_us = [bb[:8] + "_" + bb[8:] for bb in list(bc_df.values)]

            # save bb files in {SAMPLE}/bb
            pd.Series(bc_1).to_csv(output.BC_1, header=False, index=False)
            pd.Series(bc_2).to_csv(output.BC_2, header=False, index=False)
            pd.Series(bc_us).to_csv(output.BC_US, header=False, index=False)

        elif "microST" in recipes:
            bc_1 = [bb[:10] for bb in list(bc_df.values)]
            bc_1 = list(set(bc_1))

            bc_2 = [bb[10:] for bb in list(bc_df.values)]
            bc_2 = list(set(bc_2))

            bc_us = [bb[:10] + "_" + bb[10:] for bb in list(bc_df.values)]

            # save bb files in {SAMPLE}/bb
            pd.Series(bc_1).to_csv(output.BC_1, header=False, index=False)
            pd.Series(bc_2).to_csv(output.BC_2, header=False, index=False)
            pd.Series(bc_us).to_csv(output.BC_US, header=False, index=False)

        else:
            shell(
                f"""
                touch {output.BC_1}
                touch {output.BC_2}
                touch {output.BC_US}
                """
            )


# For Seeker - Insert the adapter sequence into the bead barcodes for easier barcode matching/alignment with STARsolo
rule insert_adapter_into_list:
    input:
        BC_MAP=lambda wildcards: BC_DICT[wildcards.SAMPLE],
    output:
        BC_ADAPTER="{OUTDIR}/{SAMPLE}/bb/whitelist_adapter.txt",
        BC_ADAPTER_R1="{OUTDIR}/{SAMPLE}/bb/whitelist_adapter_r1.txt",
        BC_ADAPTER_MAP="{OUTDIR}/{SAMPLE}/bb/whitelist_adapter_map.txt",
        BC_ADAPTER_R1_MAP="{OUTDIR}/{SAMPLE}/bb/whitelist_adapter_r1_map.txt",
    params:
        ADAPTER=lambda w: get_recipe_info(w, "internal.adapter"),
        R1_PRIMER=config["R1_PRIMER"],  # R1 PCR primer (Visium & Seeker)
    threads: 1
    run:
        recipes = get_recipes(wildcards, mode="concatenate")[0]

        if "seeker" in recipes:
            # load bb
            bc_df = pd.read_csv(input.BC_MAP, sep="\t", header=None)

            # split for 2 separate barcodes
            bc_1 = pd.DataFrame(bb[:8] for bb in list(bc_df[0]))
            bc_2 = pd.DataFrame(bb[8:] for bb in list(bc_df[0]))

            # Stitch bc_1, adapter, and bc_2
            bc_adapter = bc_1 + params.ADAPTER + bc_2
            bc_adapter_r1 = params.R1_PRIMER + bc_1 + params.ADAPTER + bc_2

            bc_adapter_map = pd.DataFrame(list(zip(bc_adapter[0], bc_df[1], bc_df[2])))
            bc_adapter_r1_map = pd.DataFrame(
                list(zip(bc_adapter_r1[0], bc_df[1], bc_df[2]))
            )

            # save bc files in {SAMPLE}/bb
            bc_adapter.to_csv(output.BC_ADAPTER, sep="\t", header=False, index=False)
            bc_adapter_r1.to_csv(
                output.BC_ADAPTER_R1, sep="\t", header=False, index=False
            )

            bc_adapter_map.to_csv(
                output.BC_ADAPTER_MAP, sep="\t", header=False, index=False
            )
            bc_adapter_r1_map.to_csv(
                output.BC_ADAPTER_R1_MAP, sep="\t", header=False, index=False
            )
        elif "microST" in recipes:
            # load bb
            bc_df = pd.read_csv(input.BC_MAP, sep="\t", header=None)

            # split for 2 separate barcodes
            bc_1 = pd.DataFrame(bb[0:10] for bb in list(bc_df[0]))
            bc_2 = pd.DataFrame(bb[10:] for bb in list(bc_df[0]))

            # Stitch bc_1, adapter, and bc_2
            bc_adapter = bc_1 + params.ADAPTER + bc_2
            bc_adapter_r1 = params.R1_PRIMER + bc_1 + params.ADAPTER + bc_2

            bc_adapter_map = pd.DataFrame(list(zip(bc_adapter[0], bc_df[1], bc_df[2])))
            bc_adapter_r1_map = pd.DataFrame(
                list(zip(bc_adapter_r1[0], bc_df[1], bc_df[2]))
            )

            # save bc files in {SAMPLE}/bb
            bc_adapter.to_csv(output.BC_ADAPTER, sep="\t", header=False, index=False)
            bc_adapter_r1.to_csv(
                output.BC_ADAPTER_R1, sep="\t", header=False, index=False
            )

            bc_adapter_map.to_csv(
                output.BC_ADAPTER_MAP, sep="\t", header=False, index=False
            )
            bc_adapter_r1_map.to_csv(
                output.BC_ADAPTER_R1_MAP, sep="\t", header=False, index=False
            )

        elif "decoder" in recipes:
            print("TODO")

        else:
            shell(f"touch {' '.join(output)}")
