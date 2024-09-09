# Get cell/spot/bead barcodes & UMIs
##DEPRECATED
# rule ont_umitools_extract:
#     input:
#         LOG="{OUTDIR}/{SAMPLE}/ont/misc_logs/adapter_scan_results.txt",
#         FQS=lambda w: get_fqs(w, return_type="list", mode="ONT"),
#     output:
#         R1_FQ="{OUTDIR}/{SAMPLE}/ont/barcodes_umis/{RECIPE}/umitools/umi_R1.fq.gz",
#         R2_FQ="{OUTDIR}/{SAMPLE}/ont/barcodes_umis/{RECIPE}/umitools/umi_R2.fq.gz",
#     params:
#         WHITELIST=lambda w: get_whitelist(w),
#         BC_PATTERN=lambda w: get_ont_barcode_pattern(w),
#         EXTRACT_METHOD="string",
#     log:
#         log="{OUTDIR}/{SAMPLE}/ont/barcodes_umis/{RECIPE}/umitools/extract.log",
#     resources:
#         mem="16G",
#         threads: 1,
#     run:
#         if "N" in params.BC_PATTERN:
#             shell(
#                 f"""
#                 echo "Barcode pattern: '{params.BC_PATTERN}'" > {log.log}
#                 echo "Extract method:  {params.EXTRACT_METHOD}" >> {log.log}
#                 echo "" >> {log.log}

#                 umi_tools extract \
#                     --extract-method={params.EXTRACT_METHOD} \
#                     --bc-pattern='{params.BC_PATTERN}' \
#                     --stdin={input.FQS[0]} \
#                     --read2-in={input.FQS[1]} \
#                     --stdout={output.R1_FQ} \
#                     --read2-out={output.R2_FQ} \
#                     --log2stderr \
#                 2>&1 | tee {log.log}
#                 """
#             )


# Barcode and UMI calling (custom script)
# TODO-add option for UMI-free assay
rule ont_1c_fastq_call_bc_from_adapter:
    input:
        LOG="{OUTDIR}/{SAMPLE}/ont/misc_logs/1a_adapter_scan_results.txt",
        FQS=lambda w: get_fqs(w, return_type="list", mode="ONT"),
    output:
        TSV="{OUTDIR}/{SAMPLE}/ont/barcodes_umis/{RECIPE}/read_barcodes.tsv",
    params:
        BC_PATTERN=lambda w: get_ont_barcode_pattern(w),
        BC_ADAPTERS=lambda w: get_recipe_info(w, "BC_adapter", mode="ONT"),
        BC_LENGTHS=lambda w: get_recipe_info(w, "BC_length", mode="ONT"),
        BC_OFFSETS=lambda w: get_recipe_info(w, "BC_offset", mode="ONT"),
        BC_POSITIONS=lambda w: get_recipe_info(w, "BC_position", mode="ONT"),
        BC_MISMATCHES=2,
        UMI_ADAPTERS=lambda w: get_recipe_info(w, "UMI_adapter", mode="ONT"),
        UMI_LENGTHS=lambda w: get_recipe_info(w, "UMI_length", mode="ONT"),
        UMI_OFFSETS=lambda w: get_recipe_info(w, "UMI_offset", mode="ONT"),
        UMI_POSITIONS=lambda w: get_recipe_info(w, "UMI_position", mode="ONT"),
        UMI_MISMATCHES=2,
    log:
        log="{OUTDIR}/{SAMPLE}/ont/misc_logs/{RECIPE}/1c_fastq_call_bc_from_adapter.log",
    resources:
        mem="32G",
    threads: 1
    run:
        # if any([umi_length > 0 for umi_length in params.UMI_LENGTHS]):
        shell(
            f"""
            mkdir -p $(dirname {log.log})
            python scripts/py/fastq_call_bc_umi_from_adapter.py --fq_in {input.FQS[0]} \
                --tsv_out {output.TSV} \
                --bc_adapters {params.BC_ADAPTERS} \
                --bc_lengths {params.BC_LENGTHS} \
                --bc_offsets {params.BC_OFFSETS} \
                --bc_positions {params.BC_POSITIONS} \
                --bc_mismatches {params.BC_MISMATCHES} \
                --umi_adapters {params.UMI_ADAPTERS} \
                --umi_lengths {params.UMI_LENGTHS} \
                --umi_offsets {params.UMI_OFFSETS} \
                --umi_positions {params.UMI_POSITIONS} \
                --umi_mismatches {params.UMI_MISMATCHES} \
            1> {log.log}
            """
        )
        # else:
        #     shell(
        #         f"""
        #         python scripts/py/fastq_call_bc_from_adapter.py \
        #             --fq_in {input.FQS[0]} \
        #             --tsv_out {output.TSV} \
        #             --adapters {params.BC_ADAPTERS} \
        #             --barcode_lengths {params.BC_LENGTHS} \
        #             --barcode_positions {params.BC_POSITIONS}  \
        #             --mismatches {params.BC_MISMATCHES}
        #         2>&1 | tee {log.log}
        #         """
        #     )



# Filter called read barcodes
rule ont_1c_filter_read_barcodes:
    input:
        TSV="{OUTDIR}/{SAMPLE}/ont/barcodes_umis/{RECIPE}/read_barcodes.tsv",
    output:
        TSV="{OUTDIR}/{SAMPLE}/ont/barcodes_umis/{RECIPE}/read_barcodes_filtered.tsv",
    shell:
        """
        cat {input.TSV} | grep -vP "\t-" > {output.TSV}
        """


# Correct barcodes based on white lists
# TODO- flexible input whitelist (for independent sub-barcode correction)
# TODO- recipe-specific barcode correction parameters (MAX_LEVEN, NEXT_MATCH_DIFF, CONCAT_BCS)
# TODO- option for methods which don't have UMI {params.UMI_LENGTH}
rule ont_1c_tsv_bc_correction:
    input:
        TSV="{OUTDIR}/{SAMPLE}/ont/barcodes_umis/{RECIPE}/read_barcodes_filtered.tsv",
        WHITELIST=lambda w: get_whitelist(w, return_type="string"),
        # WHITELIST="{OUTDIR}/{SAMPLE}/bc/whitelist.txt",
    output:
        TSV_SLIM="{OUTDIR}/{SAMPLE}/ont/barcodes_umis/{RECIPE}/read_barcodes_corrected.tsv",
        TSV_FULL="{OUTDIR}/{SAMPLE}/ont/barcodes_umis/{RECIPE}/read_barcodes_corrected_full.tsv",
    params:
        WHITELIST=lambda w: get_whitelist(w, return_type="string"),
        UMI_LENGTH=lambda w: get_recipe_info(w, "UMI_length", mode="ONT"),
        MAX_LEVEN=lambda w: get_recipe_info(w, "BC_max_ED", mode="ONT"),  # maximum Levenshtein distance tolerated in correction;
        NEXT_MATCH_DIFF=lambda w: get_recipe_info(w, "BC_min_ED_diff", mode="ONT"),
        K=5,  # kmer length for BC whitelist filtering; shorter value improves accuracy, extends runtime
        BC_COLUMNS=lambda w: " ".join(map(str, range(1, get_n_bcs(w) + 1))),
        CONCAT_BCS=lambda w: not get_recipe_info(w, "BC_independent", mode="ONT"),  # whether the sub-barcodes should be corrected together (SlideSeq) or separately (microST)
    log:
        log="{OUTDIR}/{SAMPLE}/ont/misc_logs/{RECIPE}/1c_tsv_bc_correction.log",
    resources:
        mem="32G",
    threads: config["CORES"]
    run:
        # if params.UMI_LENGTH > 0:
        shell(
            f"""
            python scripts/py/tsv_bc_correction_parallelized.py --tsv_in {input.TSV} \
                --tsv_out_full {output.TSV_FULL} \
                --tsv_out_slim {output.TSV_SLIM} \
                --id_column 0 \
                --bc_columns {params.BC_COLUMNS} \
                --concat_bcs {params.CONCAT_BCS} \
                --whitelist_files {input.WHITELIST} \
                --max_levens {params.MAX_LEVEN} \
                --min_next_match_diffs {params.NEXT_MATCH_DIFF} \
                --k {params.K} \
                --threads {threads} \
            1> {log.log}
            """
        )
        # else: # no UMI
        #     shell(
        #         f"""
        #         TODO
        #         """
        #     )
