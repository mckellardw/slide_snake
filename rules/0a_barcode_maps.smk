# SlideSeq-specific rules for handling the error-prone R1


# Split the bead barcodes and save whitelists
rule splitBBList:
    input:
        BB_map=lambda wildcards: BB_DICT[wildcards.SAMPLE],
    output:
        BB="{OUTDIR}/{SAMPLE}/bb/whitelist.txt",
        BB_1="{OUTDIR}/{SAMPLE}/bb/whitelist_1.txt",
        BB_2="{OUTDIR}/{SAMPLE}/bb/whitelist_2.txt",
    run:
        recipes = RECIPE_DICT[wildcards.SAMPLE]

        # get whitelist from map file
        shell(f"cut -f1 {input.BB_map} > {output.BB}")

        # split into multiple lists for Seeker
        if any(["seeker" in recipe for recipe in recipes]):
            # load bb
            bb_df = pd.read_csv(input.BB_map, sep="\t", header=None).iloc[:, 0]

            # split for 2 separate barcodes
            bb_1 = pd.DataFrame(bb[0:8] for bb in list(bb_df.values))
            bb_2 = pd.DataFrame(bb[8:14] for bb in list(bb_df.values))

            # save bb files in {SAMPLE}/bb
            bb_1.to_csv(
                output.BB_1, sep="\t", header=False, index=False
            )  # Bead barcode #1
            bb_2.to_csv(
                output.BB_2, sep="\t", header=False, index=False
            )  # Bead Barcode #2
        else:
            shell(
                f"""
                touch {output.BB_1}
                touch {output.BB_2}
                """
            )


#


# For Seeker - Insert the adapter sequence into the bead barcodes for easier barcode matching/alignment with STARsolo
rule insert_adapter_BB_list:
    input:
        BB_MAP=lambda wildcards: BB_DICT[wildcards.SAMPLE],
    output:
        BB_MAP="{OUTDIR}/{SAMPLE}/bb/whitelist_map.txt",
        BB_ADAPTER="{OUTDIR}/{SAMPLE}/bb/whitelist_adapter.txt",
        BB_ADAPTER_R1="{OUTDIR}/{SAMPLE}/bb/whitelist_adapter_r1.txt",
        BB_ADAPTER_MAP="{OUTDIR}/{SAMPLE}/bb/whitelist_adapter_map.txt",
        BB_ADAPTER_R1_MAP="{OUTDIR}/{SAMPLE}/bb/whitelist_adapter_r1_map.txt",
    params:
        ADAPTER=config["R1_INTERNAL_ADAPTER"],  # Curio R1 internal adapter
        R1=config["R1_PRIMER"],  # R1 PCR primer (Visium & Seeker)
    run:
        recipes = "".join(RECIPE_DICT[wildcards.SAMPLE])

        if "seeker" in recipes:
            # load bb
            bb_df = pd.read_csv(input.BB_MAP, sep="\t", header=None)

            # split for 2 separate barcodes
            bb_1 = pd.DataFrame(bb[0:8] for bb in list(bb_df[0]))
            bb_2 = pd.DataFrame(bb[8:14] for bb in list(bb_df[0]))

            # Stitch bb_1, adapter, and bb_2
            bb_adapter = bb_1 + params.ADAPTER + bb_2
            bb_adapter_r1 = params.R1 + bb_1 + params.ADAPTER + bb_2

            bb_adapter_map = pd.DataFrame(list(zip(bb_adapter[0], bb_df[1], bb_df[2])))
            bb_adapter_r1_map = pd.DataFrame(
                list(zip(bb_adapter_r1[0], bb_df[1], bb_df[2]))
            )

            # save bb files in {SAMPLE}/bb
            bb_adapter.to_csv(output.BB_ADAPTER, sep="\t", header=False, index=False)
            bb_adapter_r1.to_csv(
                output.BB_ADAPTER_R1, sep="\t", header=False, index=False
            )

            bb_adapter_map.to_csv(
                output.BB_ADAPTER_MAP, sep="\t", header=False, index=False
            )
            bb_adapter_r1_map.to_csv(
                output.BB_ADAPTER_R1_MAP, sep="\t", header=False, index=False
            )

        else:
            shell(f"touch {' '.join(output)}")

        shell(f"cp {input.BB_MAP} {output.BB_MAP}")
