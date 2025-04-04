# Barcode and UMI calling (custom script)
# TODO-add option for UMI-free assay
rule ont_1c_fastq_call_bc_from_adapter:
    input:
        # LOG="{OUTDIR}/{SAMPLE}/ont/logs/1a_adapter_scan_results.txt",
        FQS=lambda w: get_fqs(w, return_type="list", mode="ONT"),
    output:
        TSV="{OUTDIR}/{SAMPLE}/ont/barcodes_umis/{RECIPE}/barcodes.tsv",
        STATS_TSV="{OUTDIR}/{SAMPLE}/ont/barcodes_umis/{RECIPE}/barcode_stats.tsv",
    params:
        BC_ADAPTERS=lambda w: get_recipe_info(w, "BC_adapter", mode="ONT"),
        BC_LENGTHS=lambda w: get_recipe_info(w, "BC_length", mode="ONT"),
        BC_OFFSETS=lambda w: get_recipe_info(w, "BC_offset", mode="ONT"),
        BC_POSITIONS=lambda w: get_recipe_info(w, "BC_position", mode="ONT"),
        BC_ADAPTER_MISMATCHES=lambda w: round(
            len(get_recipe_info(w, "BC_adapter", mode="ONT")) * 0.1
        ),
        UMI_ADAPTERS=lambda w: get_recipe_info(w, "UMI_adapter", mode="ONT"),
        UMI_LENGTHS=lambda w: get_recipe_info(w, "UMI_length", mode="ONT"),
        UMI_OFFSETS=lambda w: get_recipe_info(w, "UMI_offset", mode="ONT"),
        UMI_POSITIONS=lambda w: get_recipe_info(w, "UMI_position", mode="ONT"),
        UMI_ADAPTER_MISMATCHES=lambda w: round(
            len(get_recipe_info(w, "UMI_adapter", mode="ONT")) * 0.1
        ),
    log:
        log="{OUTDIR}/{SAMPLE}/ont/barcodes_umis/{RECIPE}/1c_fastq_call_bc_from_adapter.log",
        err="{OUTDIR}/{SAMPLE}/ont/barcodes_umis/{RECIPE}/1c_fastq_call_bc_from_adapter.err",
    conda:
        f"{workflow.basedir}/envs/parasail.yml"
    resources:
        mem="32G",
    threads: 1
    shell:
        """
        mkdir -p $(dirname {log.log})
        python scripts/py/fastq_call_bc_umi_from_adapter_v2.py --fq_in {input.FQS[0]} \
            --tsv_out {output.TSV} \
            --bc_adapters {params.BC_ADAPTERS} \
            --bc_lengths {params.BC_LENGTHS} \
            --bc_offsets {params.BC_OFFSETS} \
            --bc_positions {params.BC_POSITIONS} \
            --bc_mismatches {params.BC_ADAPTER_MISMATCHES} \
            --umi_adapters {params.UMI_ADAPTERS} \
            --umi_lengths {params.UMI_LENGTHS} \
            --umi_offsets {params.UMI_OFFSETS} \
            --umi_positions {params.UMI_POSITIONS} \
            --umi_mismatches {params.UMI_ADAPTER_MISMATCHES} \
            --threads {threads} \
            --stats_out {output.STATS_TSV} \
        1> {log.log} \
        2> {log.err}
        """


# Filter called read barcodes
rule ont_1c_filter_barcodes:
    input:
        TSV="{OUTDIR}/{SAMPLE}/ont/barcodes_umis/{RECIPE}/barcodes.tsv",
    output:
        TSV="{OUTDIR}/{SAMPLE}/ont/barcodes_umis/{RECIPE}/barcodes_filtered.tsv",
    shell:
        """
        cat {input.TSV} | grep -vP "\t-" > {output.TSV}
        """


# Correct barcodes based on white lists
rule ont_1c_tsv_bc_correction:
    input:
        TSV="{OUTDIR}/{SAMPLE}/ont/barcodes_umis/{RECIPE}/barcodes_filtered.tsv",
        WHITELIST=lambda w: get_whitelist(w, return_type="list", mode="ONT"),
    output:
        TSV_SLIM="{OUTDIR}/{SAMPLE}/ont/barcodes_umis/{RECIPE}/barcodes_corrected.tsv",
        TSV_FULL="{OUTDIR}/{SAMPLE}/ont/barcodes_umis/{RECIPE}/barcodes_corrected_full.tsv",
    params:
        WHITELIST=lambda w: " ".join(get_whitelist(w, return_type="list", mode="ONT")),
        MAX_LEVEN=lambda w: get_recipe_info(w, "BC_max_ED", mode="ONT"),  # maximum Levenshtein distance tolerated in correction;
        NEXT_MATCH_DIFF=lambda w: get_recipe_info(w, "BC_min_ED_diff", mode="ONT"),
        K=5,  # kmer length for BC whitelist filtering; shorter value improves accuracy, extends runtime
        BC_COLUMNS=lambda w: " ".join(map(str, range(1, get_n_bcs(w) + 1))),
        CONCAT_BCS=lambda w: (
            "--concat_bcs" if get_recipe_info(w, "BC_concat", mode="ONT") else ""
        ),
        # whether the sub-barcodes should be corrected together (SlideSeq) or separately (microST)
    log:
        log="{OUTDIR}/{SAMPLE}/ont/barcodes_umis/{RECIPE}/1c_tsv_bc_correction.log",
        err="{OUTDIR}/{SAMPLE}/ont/barcodes_umis/{RECIPE}/1c_tsv_bc_correction.err",
    conda:
        f"{workflow.basedir}/envs/parasail.yml"
    resources:
        mem="32G",
    threads: config["CORES"]
    shell:
        """
        python scripts/py/tsv_bc_correction_parallelized.py --tsv_in {input.TSV} \
            --tsv_out_full {output.TSV_FULL} \
            --tsv_out_slim {output.TSV_SLIM} \
            --id_column 0 \
            --bc_columns {params.BC_COLUMNS} {params.CONCAT_BCS} \
            --whitelist_files {params.WHITELIST} \
            --max_levens {params.MAX_LEVEN} \
            --min_next_match_diffs {params.NEXT_MATCH_DIFF} \
            --k {params.K} \
            --threads {threads} \
        1> {log.log} \
        2> {log.err}
        """


# Summary of barcode correction results


rule ont_1c_summarize_bc_correction:
    input:
        # TSV_SLIM="{OUTDIR}/{SAMPLE}/ont/barcodes_umis/{RECIPE}/barcodes_corrected.tsv",
        TSV_FULL="{OUTDIR}/{SAMPLE}/ont/barcodes_umis/{RECIPE}/barcodes_corrected_full.tsv",
    output:
        SUMMARY="{OUTDIR}/{SAMPLE}/ont/barcodes_umis/{RECIPE}/bc_correction_stats.txt",
    params:
        N_LINES=100000,  # number of lines to summarize 
    conda:
        f"{workflow.basedir}/envs/parasail.yml"
    resources:
        mem="32G",
    threads: config["CORES"]
    shell:
        """
        python scripts/py/tsv_bc_correction_summary.py {input.TSV_FULL} > {output.SUMMARY}
        """
