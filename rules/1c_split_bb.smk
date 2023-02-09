# SPlit the bead barcodes and save whitelists
rule splitBBList:
    input:
        BB_map = lambda wildcards: BB_DICT[wildcards.sample]
    output:
        BB = "{OUTDIR}/{sample}/bb/whitelist.txt",
        BB_1 = "{OUTDIR}/{sample}/bb/whitelist_1.txt",
        BB_2 = "{OUTDIR}/{sample}/bb/whitelist_2.txt"
    run:
        #load bc
        bc_df = pd.read_csv(input.BB_map, sep="\t", header=None).iloc[:,0]

        # split for 2 separate barcodes
        bc_1 = pd.DataFrame(bc[0:8] for bc in list(bc_df.values))
        bc_2 = pd.DataFrame(bc[8:14] for bc in list(bc_df.values))

        # save bc files in {sample}/tmp
        bc_df.to_csv(output.BB, sep="\t", header=False, index=False) # Full bead barcode
        bc_1.to_csv(output.BB_1, sep="\t", header=False, index=False) # Bead barcode #1
        bc_2.to_csv(output.BB_2, sep="\t", header=False, index=False) # Bead Barcode #2
