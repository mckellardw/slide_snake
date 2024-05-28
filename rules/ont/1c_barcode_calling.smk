# Get cell/spot/bead barcodes & UMIs
rule ont_umitools_extract:
    input:
        LOG="{OUTDIR}/{SAMPLE}/ont/misc_logs/adapter_scan_results.txt",
        FQS=lambda w: get_fqs(w, return_type="list", mode="ONT"),
    output:
        R1_FQ="{OUTDIR}/{SAMPLE}/ont/barcodes_umis/{RECIPE}/umitools/umi_R1.fq.gz",
        R2_FQ="{OUTDIR}/{SAMPLE}/ont/barcodes_umis/{RECIPE}/umitools/umi_R2.fq.gz",
    params:
        WHITELIST=lambda w: get_whitelist(w),
        BC_PATTERN=lambda w: get_ont_barcode_pattern(w),
        EXTRACT_METHOD="string",
    log:
        log="{OUTDIR}/{SAMPLE}/ont/barcodes_umis/{RECIPE}/umitools/extract.log",
    threads: 1
    run:
        if "N" in params.BC_PATTERN:
            shell(
                f"""
                echo "Barcode pattern: '{params.BC_PATTERN}'" > {log.log}
                echo "Extract method:  {params.EXTRACT_METHOD}" >> {log.log}
                echo "" >> {log.log}

                {EXEC['UMITOOLS']} extract \
                    --extract-method={params.EXTRACT_METHOD} \
                    --bc-pattern='{params.BC_PATTERN}' \
                    --stdin={input.FQS[0]} \
                    --read2-in={input.FQS[1]} \
                    --stdout={output.R1_FQ} \
                    --read2-out={output.R2_FQ} \
                    --log2stderr \
                2>&1 | tee {log.log}
                """
            )


# TODO- add UMI option to script
rule ont_fastq_call_bc_from_adapter:
    input:
        LOG="{OUTDIR}/{SAMPLE}/ont/misc_logs/adapter_scan_results.txt",
        FQS=lambda w: get_fqs(w, return_type="list", mode="ONT"),
    output:
        TSV="{OUTDIR}/{SAMPLE}/ont/barcodes_umis/{RECIPE}/read_barcodes.tsv",
    params:
        BC_PATTERN=lambda w: get_ont_barcode_pattern(w),
        BC_ADAPTERS=lambda w: get_recipe_info(w, "BC.adapter", mode="ONT"),
        BC_LENGTHS=lambda w: get_recipe_info(w, "BC.length", mode="ONT"),
        BC_OFFSETS=lambda w: get_recipe_info(w, "BC.offset", mode="ONT"),
        BC_POSITIONS=lambda w: get_recipe_info(w, "BC.position", mode="ONT"),
        BC_MISMATCHES=2,
        UMI_ADAPTERS=lambda w: get_recipe_info(w, "UMI.adapter", mode="ONT"),
        UMI_LENGTHS=lambda w: get_recipe_info(w, "UMI.length", mode="ONT"),
        UMI_OFFSETS=lambda w: get_recipe_info(w, "UMI.offset", mode="ONT"),
        UMI_POSITIONS=lambda w: get_recipe_info(w, "UMI.position", mode="ONT"),
        UMI_MISMATCHES=2,
    log:
        log="{OUTDIR}/{SAMPLE}/ont/misc_logs/{RECIPE}/fastq_call_bc_from_adapter.log",
    threads: 1
    run:
        # if any([umi_length > 0 for umi_length in params.UMI_LENGTHS]):
        shell(
            f"""
            mkdir -p $(dirname {log.log})
            python scripts/py/fastq_call_bc_umi_from_adapter.py \
                --fq_in {input.FQS[0]} \
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
            2>&1 | tee {log.log}
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



# Correct barcodes based on white lists
rule ont_tsv_bc_correction:
    input:
        TSV="{OUTDIR}/{SAMPLE}/ont/barcodes_umis/{RECIPE}/read_barcodes.tsv",
        # WHITELIST=lambda w: get_whitelist(w, return_type="list"),
        WHITELIST="{OUTDIR}/{SAMPLE}/bc/whitelist.txt",
    output:
        TSV_SLIM="{OUTDIR}/{SAMPLE}/ont/barcodes_umis/{RECIPE}/read_barcodes_corrected.tsv",
        TSV_FULL="{OUTDIR}/{SAMPLE}/ont/barcodes_umis/{RECIPE}/read_barcodes_corrected_full.tsv",
    params:
        # BC_PATTERN=lambda w: get_ont_barcode_pattern(w),
        MAX_HAM=2,
        BC_COLUMNS=lambda w: " ".join(map(str, range(1, get_n_bcs(w) + 1))),
    log:
        log="{OUTDIR}/{SAMPLE}/ont/misc_logs/{RECIPE}/tsv_bc_correction.log",
    threads: 1
    run:
        # if any([umi_length > 0 for umi_length in params.UMI_LENGTHS]):
        shell(
            f"""
            python scripts/py/tsv_bc_correction.py \
                --tsv_in {input.TSV} \
                --tsv_out_full {output.TSV_FULL} \
                --tsv_out_slim {output.TSV_SLIM} \
                --id_column 0 \
                --bc_columns {params.BC_COLUMNS} \
                --concat_bcs True \
                --whitelist_files {input.WHITELIST} \
                --max_ham {params.MAX_HAM} \
            2>&1 | tee {log.log}
            """
        )
        # else: # no UMI
        #     shell(
        #         f"""
        #         TODO
        #         """
        #     )
