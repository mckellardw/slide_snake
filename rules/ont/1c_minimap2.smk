
rule ont_align_minimap2:
    input:
        FQ = "{OUTDIR}/{SAMPLE}/tmp/ont/cut_R2.fq.gz",
    output:
        SAM_TMP=temp("{OUTDIR}/{SAMPLE}/ont/minimap2/tmp.sam"),
    params:
        ref = config["REF_GENOME_FASTA"],
        chrom_sizes = config["REF_CHROM_SIZES"],
        bed = config["REF_GENES_BED"],
        flags = config["RESOURCES_MM2_FLAGS"],
    log:
        log = "{OUTDIR}/{SAMPLE}/ont/minimap2/minimap2.log"
    threads: 
        config["CORES"]
    # resources:
    #     mem_gb=get_split_ont_align_mem_gb,
    # conda:
    #     "../envs/minimap2.yml"
    run:
        shell(
            f"""
            {EXEC['MINIMAP2']} \
                -ax splice \
                -uf \
                --MD \
                -t {threads} \
                --junc-bed {params.bed} \
                --secondary=no \
                {params.flags} \
                {params.ref} \
                {input.FQ} \
            2> {log.log} \
            > {output.SAM_TMP}
            """
        )
#


rule ont_sort_index_output:
    input:
        SAM="{OUTDIR}/{SAMPLE}/ont/minimap2/tmp.sam"
    output:
        # BAM_UNSORT_TMP=temp("{OUTDIR}/{SAMPLE}/ont/tmp_unsort.sam"),
        BAM="{OUTDIR}/{SAMPLE}/ont/minimap2/sorted.bam",
        BAI="{OUTDIR}/{SAMPLE}/ont/minimap2/sorted.bam.bai",
    params:
        ref = config["REF_GENOME_FASTA"]
    # threads: 
    #     config["CORES"]
    run:
        shell(
            f"""
            {EXEC['SAMTOOLS']} sort \
                --reference {params.ref} \
                -O BAM \
                -o {output.BAM} \
                {input.SAM}
            
            {EXEC['SAMTOOLS']} index {output.BAM}
            """
        )
#


rule ont_extract_barcodes:
    input:
        BAM = "{OUTDIR}/{SAMPLE}/ont/minimap2/sorted.bam",
        BB_WHITELIST = "{OUTDIR}/{SAMPLE}/bb/whitelist.txt"
    output:
        BAM = "{OUTDIR}/{SAMPLE}/ont/minimap2/sorted.bam", #temp()
        TSV = "{OUTDIR}/{SAMPLE}/ont/minimap2/bc_counts.tsv",
    params:
        KIT = "3prime", #['3prime', '5prime', 'multiome']
        adapter1_suff_length = config["BARCODE_ADAPTER1_SUFF_LENGTH"],
        barcode_min_qv = config["BARCODE_MIN_QUALITY"],
        barcode_length = lambda w: get_barcode_length(w),
        umi_length = lambda w: get_umi_length(w),
    threads: 
        config["CORES"]
    # conda:
    #     "../envs/barcodes.yml"
    run:
        shell(
            f"""
            python scripts/py/extract_barcode.py \
                -t {threads} \
                --kit {params.kit} \
                --adapter1_suff_length {params.adapter1_suff_length} \
                --min_barcode_qv {params.barcode_min_qv} \
                --barcode_length {params.barcode_length} \
                --umi_length {params.umi_length} \
                --output_bam {output.BAM} \
                --output_barcodes {output.TSV} \
                {input.BAM} {input.BB_WHITELIST}
            """
        )
#

# rule cleanup_headers_1:
#     input:
#         BAM_BC_UNCORR_TMP,
#     output:
#         bam=temp(BAM_BC_UNCORR),
#         bai=temp(BAM_BC_UNCORR_BAI),
#     # conda:
#     #     "../envs/{EXEC['SAMTOOLS']}.yml"
#     run:
#         shell(
#             f"""
#             {EXEC['SAMTOOLS']} reheader \
#                 --no-PG \
#                 -c 'grep -v ^@PG' \
#                 {input} \
#             > {output.bam}
            
#             {EXEC['SAMTOOLS']} index {output.bam}
#             """
#         )


# rule assign_barcodes:
#     input:
#         bam=CHROM_BAM_BC_UNCORR,
#         bai=CHROM_BAM_BC_UNCORR_BAI,
#         whitelist=BARCODE_WHITELIST,
#     output:
#         bam=temp(CHROM_BAM_BC_TMP),
#         counts=CHROM_ASSIGNED_BARCODE_COUNTS,
#     params:
#         max_ed=config["BARCODE_MAX_ED"],
#         min_ed_diff=config["BARCODE_MIN_ED_DIFF"],
#         # kit=lambda w: sample_sheet.loc[w.run_id, "kit_name"],
#         KIT = "3prime", #['3prime', '5prime', 'multiome']
#         adapter1_suff_length=config["BARCODE_ADAPTER1_SUFF_LENGTH"],
#         barcode_length=lambda w: get_barcode_length(w),
#         umi_length=lambda w: get_umi_length(w),
#         # barcode_length=config["READ_STRUCTURE_BARCODE_LENGTH"],
#         # umi_length=config["READ_STRUCTURE_UMI_LENGTH"],
#     threads: 1
#     # conda:
#     #     "../envs/barcodes.yml"
#     run:
#         shell(
#             f"""
#             touch {input.bai}

#             python scripts/py/assign_barcodes.py \
#                 -t {threads} \
#                 --output_bam {output.bam} \
#                 --output_counts {output.counts} \
#                 --max_ed {params.max_ed} \
#                 --min_ed_diff {params.min_ed_diff} \
#                 --kit {params.kit} \
#                 --adapter1_suff_length {params.adapter1_suff_length} \
#                 --barcode_length {params.barcode_length} \
#                 --umi_length {params.umi_length} \
#                 {input.bam} {input.whitelist}
#             """
#         )


# Generate count matrix w/ umi-tools for miRNAs
# rule umitools_count_ONT:
#     input:
#         BAM = '{OUTDIR}/{SAMPLE}/ont/aligned.sorted.tagged.bam',
#         BAI = '{OUTDIR}/{SAMPLE}/ont/aligned.sorted.tagged.bam.bai'
#     output:        
#         COUNTS = '{OUTDIR}/{SAMPLE}/ont/counts.tsv.gz'
#     params:
#         OUTDIR = config['OUTDIR']
#     log:
#         log = '{OUTDIR}/{SAMPLE}/ont/count.log'
#     threads:
#         1
#     run:
#         shell(
#             f"""
#             {EXEC['UMITOOLS']} count \
#                 --extract-umi-method=tag \
#                 --per-gene \
#                 --per-cell \
#                 --cell-tag=CB \
#                 --gene-tag=GN \
#                 --umi-tag=UB \
#                 --log={log.log} \
#                 -I {input.BAM} \
#                 -S {output.COUNTS}
#             """
#         )