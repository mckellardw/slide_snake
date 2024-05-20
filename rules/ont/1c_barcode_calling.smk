# Get cell/spot/bead barcodes & UMIs
rule ont_umitools_extract:
    input:
        LOG="{OUTDIR}/{SAMPLE}/ont/misc_logs/adapter_scan_results.txt",
        FQS=lambda w: get_fqs(w, return_type="list", mode="ONT"),
    output:
        R1_FQ="{OUTDIR}/{SAMPLE}/ont/umitools/{RECIPE}/umi_R1.fq.gz",
        R2_FQ="{OUTDIR}/{SAMPLE}/ont/umitools/{RECIPE}/umi_R2.fq.gz",
    params:
        WHITELIST=lambda w: get_whitelist(w),
        BC_PATTERN=lambda w: get_ont_barcode_pattern(w),
        EXTRACT_METHOD="string",
    log:
        log="{OUTDIR}/{SAMPLE}/ont/umitools/{RECIPE}/extract.log",
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

#TODO
rule ont_fastq_call_bc_from_adapter:
    input:
        LOG="{OUTDIR}/{SAMPLE}/ont/misc_logs/adapter_scan_results.txt",
        FQS=lambda w: get_fqs(w, return_type="list", mode="ONT"),
    output:
        TSV="{OUTDIR}/{SAMPLE}/ont/{RECIPE}/read_barcodes.tsv"
    params:
        WHITELIST=lambda w: get_whitelist(w),
        BC_PATTERN=lambda w: get_ont_barcode_pattern(w),
        EXTRACT_METHOD="string",
    log:
        log="{OUTDIR}/{SAMPLE}/ont/misc_logs/{RECIPE}/fastq_call_bc_from_adapter.log",
    threads: 1
    run:
        if "N" in params.BC_PATTERN:
            shell(
                f"""
                python scripts/py/fastq_call_bc_from_adapter.py \
                    --fq_in sandbox/lig3_R1.fq \
                    --tsv_out sandbox/barcodes.tsv \
                    --adapters TCCACGTGCTTGAG TCCACGTGCTTGAG \
                    --barcode_lengths 10 10 \
                    --barcode_positions left right \
                    --mismatches 1 1
                2>&1 | tee {log.log}
                """
            )
