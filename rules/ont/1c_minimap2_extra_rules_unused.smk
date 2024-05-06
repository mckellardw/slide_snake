

# Old ONT code - barcode is left in the read during alignment which seems dumb?
# rule ont_extract_barcodes_from_R1:
#     input:
#         # FQ = "{OUTDIR}/{SAMPLE}/tmp/ont/cut_R1.fq.gz",
#         BAM = "{OUTDIR}/{SAMPLE}/ont/minimap2/sorted.bam",
#         BB_WHITELIST = "{OUTDIR}/{SAMPLE}/bb/whitelist.txt",
#         BB_1 = "{OUTDIR}/{SAMPLE}/bb/whitelist_1.txt",
#         BB_2 = "{OUTDIR}/{SAMPLE}/bb/whitelist_2.txt",
#         BB_ADAPTER = "{OUTDIR}/{SAMPLE}/bb/whitelist_adapter.txt",
#         BB_ADAPTER_R1 = "{OUTDIR}/{SAMPLE}/bb/whitelist_adapter_r1.txt",
#     output:
#         BAM = "{OUTDIR}/{SAMPLE}/ont/minimap2/{RECIPE}/sorted_bc.bam", #temp()
#         TSV = "{OUTDIR}/{SAMPLE}/ont/minimap2/{RECIPE}/bc_counts.tsv",
#     params:
#         KIT = "3prime", #['3prime', '5prime']
#         adapter1_suff_length = config["BARCODE_ADAPTER1_SUFF_LENGTH"],
#         barcode_min_qv = config["BARCODE_MIN_QUALITY"],
#         WHITELIST=lambda w: get_whitelist(w),
#         # barcode_length = lambda w: get_barcode_length(w),
#         umi_length = lambda w: get_umi_length(w),
#     threads:
#         config["CORES"]
#     log:
#         log = "{OUTDIR}/{SAMPLE}/ont/minimap2/{RECIPE}/extract_barcodes.log",
#     # conda:
#     #     "../envs/barcodes.yml"
#     run:
#         # recipe = RECIPE_ONT_DICT[wildcards.SAMPLE]
#         recipe = wildcards.RECIPE
#         shell(
#             f"""
#             echo "UMI length: {params.umi_length}" >> {log.log}

#             python scripts/py/extract_barcode.py \
#                 -t {threads} \
#                 --kit {params.KIT} \
#                 --adapter1_suff_length {params.adapter1_suff_length} \
#                 --min_barcode_qv {params.barcode_min_qv} \
#                 --umi_length {params.umi_length} \
#                 --output_bam {output.BAM} \
#                 --output_barcodes {output.TSV} \
#                 {input.BAM} {params.WHITELIST} \
#             2>> {log.log}
#             """
#         )
# barcode_length = get_barcode_length(wildcards)
# barcode_length = len(open(whitelist).readline())
# echo "Barcode length: {barcode_length}" > {log.log}
# --barcode_length {barcode_length} \



# rule ont_assign_genes:
#     input:
#         bed=CHROM_BED_BC,
#         chrom_gtf=CHROM_GTF,
#     output:
#         GENES=CHROM_TSV_GENE_ASSIGNS,
#     params:
#         minQV=config["GENE_ASSIGNS_MINQV"],
#     threads:
#     1
#     conda:
#         "../envs/assign_genes.yml"
#     shell:
#         f"""
#         python {SCRIPT_DIR}/assign_genes.py \
#             --output {output.GENES} \
#             {input.bed} {input.chrom_gtf}
#         """

# rule ont_add_gene_tags_to_bam:
#     input:
#         BAM="{OUTDIR}/{SAMPLE}/ont/minimap2/{RECIPE}/sorted_bc_corrected.bam",
#         BAI="{OUTDIR}/{SAMPLE}/ont/minimap2/{RECIPE}/sorted_bc_corrected.bam.bai",
#         GENES=CHROM_TSV_GENE_ASSIGNS,
#     output:
#         bam=temp(CHROM_BAM_BC_GENE_TMP),
#     conda:
#         "../envs/umis.yml"
#     shell:
#         "touch {input.bai}; "
#         "python {SCRIPT_DIR}/add_gene_tags.py "
#         "--output {output.bam} "
#         "{input.bam} {input.genes}"
# #

# featureCounts details:
## -s 1 = stranded
## -L = long-read mode
## -f = feature-level labelling (as in, not 'gene'-level, but 'transcript' level)
# rule ont_featureCounts:
#     input:
#         BAM="{OUTDIR}/{SAMPLE}/ont/minimap2/{RECIPE}/sorted_bc.bam",
#         BAI="{OUTDIR}/{SAMPLE}/ont/minimap2/{RECIPE}/sorted_bc.bam.bai",
#     output:
#         BAM=temp("{OUTDIR}/{SAMPLE}/ont/minimap2/{RECIPE}/sorted_bc.bam.featureCounts.bam"),
#         FEAT="{OUTDIR}/{SAMPLE}/ont/minimap2/{RECIPE}/featureCounts.tsv",
#     params:
#         GTF = lambda wildcards: GTF_DICT[wildcards.SAMPLE],
#         EXTRA_FLAGS = lambda wildcards: RECIPE_SHEET["featureCounts.extra"][wildcards.RECIPE],
#         MIN_TEMPLATE_LENGTH = 10,
#         MAX_TEMPLATE_LENGTH = 10000
#     log:
#         log = "{OUTDIR}/{SAMPLE}/ont/minimap2/{RECIPE}/featureCounts.log"
#     threads:
#         1 # long reads can only run single-threaded
#     # conda:
#     run:
#         shell(
#             f"""
#             {EXEC['FEATURECOUNTS']} \
#                 -a {params.GTF} \
#                 -o {output.FEAT} \
#                 -L \
#                 -s 1 \
#                 -f \
#                 -d {params.MIN_TEMPLATE_LENGTH} \
#                 -D {params.MAX_TEMPLATE_LENGTH} \
#                 -t 'transcript' \
#                 -g 'transcript_id' \
#                 -T {threads} \
#                 --donotsort \
#                 -R BAM {params.EXTRA_FLAGS} \
#                 {input.BAM} \
#             2> {log.log}
#             """
#         )
