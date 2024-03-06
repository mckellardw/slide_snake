
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
        BAM=temp("{OUTDIR}/{SAMPLE}/ont/minimap2/sorted.bam"),
        BAI=temp("{OUTDIR}/{SAMPLE}/ont/minimap2/sorted.bam.bai"),
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


rule ont_extract_barcodes_from_R1:
    input:
        # FQ = "{OUTDIR}/{SAMPLE}/tmp/ont/cut_R1.fq.gz",
        BAM = "{OUTDIR}/{SAMPLE}/ont/minimap2/sorted.bam",
        BB_WHITELIST = "{OUTDIR}/{SAMPLE}/bb/whitelist.txt",
        BB_1 = "{OUTDIR}/{SAMPLE}/bb/whitelist_1.txt",
        BB_2 = "{OUTDIR}/{SAMPLE}/bb/whitelist_2.txt",
        BB_ADAPTER = "{OUTDIR}/{SAMPLE}/bb/whitelist_adapter.txt"
    output:
        BAM = "{OUTDIR}/{SAMPLE}/ont/minimap2/sorted_bc.bam", #temp()
        TSV = "{OUTDIR}/{SAMPLE}/ont/minimap2/bc_counts.tsv",
    params:
        KIT = "3prime", #['3prime', '5prime', 'multiome']
        adapter1_suff_length = config["BARCODE_ADAPTER1_SUFF_LENGTH"],
        barcode_min_qv = config["BARCODE_MIN_QUALITY"],
        # barcode_length = lambda w: get_barcode_length(w),
        # umi_length = lambda w: get_umi_length(w),
    threads: 
        config["CORES"]
    # conda:
    #     "../envs/barcodes.yml"
    run:
        recipe = RECIPE_DICT[wildcards.SAMPLE]
        # recipe = wildcards.RECIPE
        
        #param handling for different SlideSeq R1 strategies
        if "stomics" in recipe:
            whitelist = input.BB_WHITELIST
        elif "noTrim" in recipe or "matchLinker" in recipe:
            whitelist = f"{input.BB_1} {input.BB_2}"
        elif "internalTrim" in recipe:
            whitelist = input.BB_WHITELIST
        elif "adapterInsert" in recipe:
            whitelist = input.BB_ADAPTER
        else:
            whitelist = input.BB_WHITELIST

        # barcode_length = get_barcode_length(wildcards)
        barcode_length = len(open(whitelist).readline())
        umi_length = get_umi_length(wildcards)

        shell(
            f"""
            python scripts/py/extract_barcode.py \
                -t {threads} \
                --kit {params.KIT} \
                --adapter1_suff_length {params.adapter1_suff_length} \
                --min_barcode_qv {params.barcode_min_qv} \
                --barcode_length {barcode_length} \
                --umi_length {umi_length} \
                --output_bam {output.BAM} \
                --output_barcodes {output.TSV} \
                {input.BAM} {whitelist}
            """
        )
#

rule ont_index_bc_output:
    input:
        BAM="{OUTDIR}/{SAMPLE}/ont/minimap2/sorted_bc.bam"
    output:
        BAI="{OUTDIR}/{SAMPLE}/ont/minimap2/sorted_bc.bam.bai",
    params:
        ref = config["REF_GENOME_FASTA"]
    threads:
        1
        # config["CORES"]
    run:
        shell(
            f"""
            {EXEC['SAMTOOLS']} index {input.BAM}
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

#TODO - need to fix the script so this can be parallelized
rule assign_barcodes:
    input:
        BAM = "{OUTDIR}/{SAMPLE}/ont/minimap2/sorted_bc.bam",
        BAI = "{OUTDIR}/{SAMPLE}/ont/minimap2/sorted_bc.bam.bai",
        BB_WHITELIST = "{OUTDIR}/{SAMPLE}/bb/whitelist.txt",
        BB_1 = "{OUTDIR}/{SAMPLE}/bb/whitelist_1.txt",
        BB_2 = "{OUTDIR}/{SAMPLE}/bb/whitelist_2.txt",
        BB_ADAPTER = "{OUTDIR}/{SAMPLE}/bb/whitelist_adapter.txt",
    output:
        BAM = "{OUTDIR}/{SAMPLE}/ont/minimap2/sorted_bc_corrected.bam",
        COUNTS="{OUTDIR}/{SAMPLE}/ont/minimap2/counts.tsv",
    params:
        max_ed=config["BARCODE_MAX_ED"],
        min_ed_diff=config["BARCODE_MIN_ED_DIFF"],
        # kit=lambda w: sample_sheet.loc[w.run_id, "kit_name"],
        KIT = "3prime", #['3prime', '5prime', 'multiome']
        adapter1_suff_length=config["BARCODE_ADAPTER1_SUFF_LENGTH"],
    threads: 
        # config["CORES"] #TODO
        1
    # conda:
    #     "../envs/barcodes.yml"
    run:
        recipe = RECIPE_DICT[wildcards.SAMPLE]
        # recipe = wildcards.RECIPE

        barcode_length = get_barcode_length(wildcards)
        umi_length = get_umi_length(wildcards)
        
        #param handling for different SlideSeq R1 strategies
        if "stomics" in recipe:
            whitelist = input.BB_WHITELIST
        elif "noTrim" in recipe or "matchLinker" in recipe:
            whitelist = f"{input.BB_1} {input.BB_2}"
        elif "internalTrim" in recipe:
            whitelist = input.BB_WHITELIST
        elif "adapterInsert" in recipe:
            whitelist = input.BB_ADAPTER
        else:
            whitelist = input.BB_WHITELIST

        shell(
            f"""
            python scripts/py/assign_barcodes_noChromSplit.py \
                -t {threads} \
                --output_bam {output.BAM} \
                --output_counts {output.COUNTS} \
                --max_ed {params.max_ed} \
                --min_ed_diff {params.min_ed_diff} \
                --kit {params.KIT} \
                --adapter1_suff_length {params.adapter1_suff_length} \
                --barcode_length {barcode_length} \
                --umi_length {umi_length} \
                {input.BAM} {whitelist}
            """
        )


# Generate count matrix w/ umi-tools 
rule ont_umitools_count:
    input:
        BAM = "{OUTDIR}/{SAMPLE}/ont/minimap2/sorted_bc.bam",
        BAI = "{OUTDIR}/{SAMPLE}/ont/minimap2/sorted_bc.bam.bai",
    output:        
        COUNTS = '{OUTDIR}/{SAMPLE}/ont/minimap2/umitools_counts.tsv.gz'
    params:
        CELL_TAG="CB",
        GENE_TAG="GN",
        UMI_TAG="UB"
    log:
        log = '{OUTDIR}/{SAMPLE}/ont/minimap2/count.log'
    threads:
        1
    run:
        shell(
            f"""
            {EXEC['UMITOOLS']} count \
                --extract-umi-method=tag \
                --per-gene \
                --per-cell \
                --cell-tag={params.CELL_TAG} \
                --gene-tag={params.GENE_TAG}  \
                --umi-tag={params.UMI_TAG}  \
                --log={log.log} \
                -I {input.BAM} \
                -S {output.COUNTS} 
            """
        )