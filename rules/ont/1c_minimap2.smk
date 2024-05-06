# Get cell/spot/bead barcodes & UMIs
rule ont_umitools_extract:
    input:
        # R1_FQ="{OUTDIR}/{SAMPLE}/tmp/ont/cut_R1.fq.gz",
        R1_FQ="{OUTDIR}/{SAMPLE}/tmp/ont/cut_trim_R1.fq.gz",
        R2_FQ="{OUTDIR}/{SAMPLE}/tmp/ont/cut_R2.fq.gz",
        BB_WHITELIST="{OUTDIR}/{SAMPLE}/bb/whitelist.txt",
        BB_1="{OUTDIR}/{SAMPLE}/bb/whitelist_1.txt",
        BB_2="{OUTDIR}/{SAMPLE}/bb/whitelist_2.txt",
        BB_ADAPTER="{OUTDIR}/{SAMPLE}/bb/whitelist_adapter.txt",
        BB_ADAPTER_R1="{OUTDIR}/{SAMPLE}/bb/whitelist_adapter_r1.txt",
    output:
        R1_FQ="{OUTDIR}/{SAMPLE}/ont/umitools/{RECIPE}/umi_R1.fq.gz",
        R2_FQ="{OUTDIR}/{SAMPLE}/ont/umitools/{RECIPE}/umi_R2.fq.gz",
    params:
        WHITELIST=lambda w: get_whitelist(w),
        BC_PATTERN=lambda w: get_ont_barcode_pattern(w),
        # BC_PATTERN="CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCNNNNNNN" #TODO
        EXTRACT_METHOD="string",
    log:
        log="{OUTDIR}/{SAMPLE}/ont/umitools/{RECIPE}/extract.log",
    threads: 1
    run:
        shell(
            f"""
            echo "Barcode pattern: '{params.BC_PATTERN}'" > {log.log}
            echo "Extract method:  {params.EXTRACT_METHOD}" >> {log.log}
            echo "" >> {log.log}

            {EXEC['UMITOOLS']} extract \
                --extract-method={params.EXTRACT_METHOD} \
                --bc-pattern='{params.BC_PATTERN}' \
                --stdin={input.R1_FQ} \
                --read2-in={input.R2_FQ} \
                --stdout={output.R1_FQ} \
                --read2-out={output.R2_FQ} \
                --log2stderr \
            2>> {log.log}
            """
        )
        # --error-correct-cell \
        # --whitelist={whitelist} \
        # --error-correct-cell \



# Align w/ minimap2
## minimap2 docs - https://lh3.github.io/minimap2/minimap2.html
rule ont_align_minimap2_genome:
    input:
        FQ="{OUTDIR}/{SAMPLE}/ont/umitools/{RECIPE}/umi_R2.fq.gz",
    output:
        SAM_TMP=temp("{OUTDIR}/{SAMPLE}/ont/minimap2/{RECIPE}/tmp.sam"),
    params:
        EXTRA_FLAGS=lambda wildcards: RECIPE_SHEET["minimap2.extra"][wildcards.RECIPE],
        ref=config["REF_GENOME_FASTA"],
        chrom_sizes=config["REF_CHROM_SIZES"],
        bed=config["REF_GENES_BED"],
        flags=config["RESOURCES_MM2_FLAGS"],
    log:
        log="{OUTDIR}/{SAMPLE}/ont/minimap2/{RECIPE}/minimap2.log",
    threads: config["CORES"]
    # resources:
    #     mem_gb=get_split_ont_align_mem_gb,
    # conda:
    #     f"{workflow.basedir}/envs/minimap2.yml" #TODO
    run:
        shell(
            f"""
            mkdir -p $(dirname {output.SAM_TMP})

            echo "Extra flags: {params.EXTRA_FLAGS}" > {log.log} 
            echo "" >> {log.log} 

            {EXEC['MINIMAP2']} \
                -ax splice \
                -uf \
                --MD \
                -t {threads} \
                --junc-bed {params.bed} \
                {params.flags} {params.EXTRA_FLAGS} \
                {params.ref} \
                {input.FQ} \
            2>> {log.log} \
            > {output.SAM_TMP}
            """
        )


#
# if "total" in wildcards.RECIPE:
#     EXTRA_FLAGS = "-M --secondary=yes"
# else:
#     EXTRA_FLAGS = "--secondary=no"


rule ont_sort_index_output:
    input:
        SAM="{OUTDIR}/{SAMPLE}/ont/minimap2/{RECIPE}/tmp.sam",
    output:
        # BAM_UNSORT_TMP=temp("{OUTDIR}/{SAMPLE}/ont/tmp_unsort.sam"),
        BAM=temp("{OUTDIR}/{SAMPLE}/ont/minimap2/{RECIPE}/sorted.bam"),
        # BAI=temp("{OUTDIR}/{SAMPLE}/ont/minimap2/{RECIPE}/sorted.bam.bai"),
    params:
        ref=config["REF_GENOME_FASTA"],
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
            
            """
        )
            # {EXEC['SAMTOOLS']} index {output.BAM}



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


# Get UMI & CB from R1 fastq
rule ont_extract_barcodes_from_R1:
    input:
        BAM="{OUTDIR}/{SAMPLE}/ont/minimap2/{RECIPE}/sorted.bam",
    output:
        BAM="{OUTDIR}/{SAMPLE}/ont/minimap2/{RECIPE}/sorted_bc.bam",
    params:
        BARCODE_TAG="CR",  # uncorrected!
        UMI_TAG="UR",  # uncorrected!
    threads: 1
    log:
        log="{OUTDIR}/{SAMPLE}/ont/minimap2/{RECIPE}/extract_barcodes.log",
    run:
        shell(
            f"""
            python scripts/py/readID2tags.py \
                --input {input.BAM} \
                --output {output.BAM} \
                -c {params.BARCODE_TAG} \
                -y {params.UMI_TAG} \
            2>> {log.log}
            """
        )


# rule ont_index_bc_output:
#     input:
#         BAM="{OUTDIR}/{SAMPLE}/ont/minimap2/{RECIPE}/sorted_bc.bam",
#     output:
#         BAI="{OUTDIR}/{SAMPLE}/ont/minimap2/{RECIPE}/sorted_bc.bam.bai",
#     threads: 1
#     run:
#         shell(
#             f"""
#             {EXEC['SAMTOOLS']} index {input.BAM}
#             """
#         )


# Note - run time increases (substantially) with BC length
rule ont_error_correct_barcodes:
    input:
        BAM="{OUTDIR}/{SAMPLE}/ont/minimap2/{RECIPE}/sorted_bc.bam",
        BAI="{OUTDIR}/{SAMPLE}/ont/minimap2/{RECIPE}/sorted_bc.bam.bai",
        BB_WHITELIST="{OUTDIR}/{SAMPLE}/bb/whitelist.txt",
        BB_1="{OUTDIR}/{SAMPLE}/bb/whitelist_1.txt",
        BB_2="{OUTDIR}/{SAMPLE}/bb/whitelist_2.txt",
        BB_ADAPTER="{OUTDIR}/{SAMPLE}/bb/whitelist_adapter.txt",
        BB_ADAPTER_R1="{OUTDIR}/{SAMPLE}/bb/whitelist_adapter_r1.txt",
    output:
        BAM="{OUTDIR}/{SAMPLE}/ont/minimap2/{RECIPE}/sorted_bc_corrected.bam",
        COUNTS="{OUTDIR}/{SAMPLE}/ont/minimap2/{RECIPE}/counts.tsv",
    params:
        max_ed=config["BARCODE_MAX_ED"],
        min_ed_diff=config["BARCODE_MIN_ED_DIFF"],
        adapter1_suff_length=config["BARCODE_ADAPTER1_SUFF_LENGTH"],
        # kit=lambda w: sample_sheet.loc[w.run_id, "kit_name"],
        KIT="3prime",  #['3prime', '5prime', 'multiome']
        # UMI_LENGTH=lambda w: get_umi_length(w),
        UMI_LENGTH=lambda w: get_recipe_info(w.RECIPE, "UMI.length"),
        WHITELIST=lambda w: get_whitelist(w),
    threads: config["CORES"]
    # 1
    log:
        log="{OUTDIR}/{SAMPLE}/ont/minimap2/{RECIPE}/correct_barcodes.log",
    # conda:
    #     "../envs/barcodes.yml"
    run:
        shell(
            f"""
            echo "UMI length:        {params.umi_length}" > {log.log}
            echo "Barcode whitelist: {params.WHITELIST}" >> {log.log}

            python scripts/py/barcodes_correct_parallelized.py \
                -t {threads} \
                --output_bam {output.BAM} \
                --output_counts {output.COUNTS} \
                --max_ed {params.max_ed} \
                --min_ed_diff {params.min_ed_diff} \
                --kit {params.KIT} \
                --adapter1_suff_length {params.adapter1_suff_length} \
                --umi_length {params.umi_length} \
                {input.BAM} {params.WHITELIST} \
            2>> {log.log}
            """
        )
        # barcode_length = get_barcode_length(wildcards)
        # barcode_length = len(open(whitelist).readline())
        # umi_length = get_umi_length(wildcards)
        # echo "Barcode length: {barcode_length}" > {log.log}
        #     --barcode_length {barcode_length} \



# rule ont_index_bc_corrected_output:
#     input:
#         BAM="{OUTDIR}/{SAMPLE}/ont/minimap2/{RECIPE}/sorted_bc_corrected.bam",
#     output:
#         BAI="{OUTDIR}/{SAMPLE}/ont/minimap2/{RECIPE}/sorted_bc_corrected.bam.bai",
#     threads: 1
#     # config["CORES"]
#     run:
#         shell(
#             f"""
#             {EXEC['SAMTOOLS']} index {input.BAM}
#             """
#         )


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


rule ont_featureCounts:
    input:
        BAM="{OUTDIR}/{SAMPLE}/ont/minimap2/{RECIPE}/sorted_bc_corrected.bam",
        BAI="{OUTDIR}/{SAMPLE}/ont/minimap2/{RECIPE}/sorted_bc_corrected.bam.bai",
    output:
        TSV="{OUTDIR}/{SAMPLE}/ont/minimap2/{RECIPE}/sorted_bc_corrected.bam.featureCounts",
        FEAT="{OUTDIR}/{SAMPLE}/ont/minimap2/{RECIPE}/featureCounts.tsv",
    params:
        GTF=lambda wildcards: GTF_DICT[wildcards.SAMPLE],
        EXTRA_FLAGS=lambda wildcards: RECIPE_SHEET["featureCounts.extra"][
            wildcards.RECIPE
        ],
        MIN_TEMPLATE_LENGTH=10,
        MAX_TEMPLATE_LENGTH=10000,
    log:
        log="{OUTDIR}/{SAMPLE}/ont/minimap2/{RECIPE}/featureCounts.log",
    threads: 1  # long reads can only run single-threaded
    # conda:
    run:
        shell(
            f"""
            {EXEC['FEATURECOUNTS']} \
                -a {params.GTF} \
                -o {output.FEAT} \
                -L \
                -s 1 \
                -f \
                -d {params.MIN_TEMPLATE_LENGTH} \
                -D {params.MAX_TEMPLATE_LENGTH} \
                -t 'transcript' \
                -g 'transcript_id' \
                -T {threads} \
                -R CORE {params.EXTRA_FLAGS} \
                {input.BAM} \
            2> {log.log}
            """
        )


#
# --donotsort \
# −−sortReadsByCoordinates \


# rule salmon_quant:
#     input:
#         BAM="{OUTDIR}/{SAMPLE}/ont/minimap2/{RECIPE}/sorted_bc.bam",
#         BAI="{OUTDIR}/{SAMPLE}/ont/minimap2/{RECIPE}/sorted_bc.bam.bai",
#         INDEX="{PATH_TO_SALMON_INDEX}"
#     output:
#         QUANT="{OUTDIR}/{SAMPLE}/ont/minimap2/{RECIPE}/quant.sf"
#     params:
#         LIBTYPE = "A", # Automatic library type detection
#         VALIDATE_MAPPINGS = True, # Validate mappings
#         GC_BIAS = True, # Correct GC bias
#         NUM_GIBBS_SAMPLES = 20 # Number of Gibbs samples
#     log:
#         log = "{OUTDIR}/{SAMPLE}/ont/minimap2/{RECIPE}/salmon_quant.log"
#     threads:
#         config["CORES"]
#     conda:
#          f"{workflow.basedir}/envs/salmon.yml"
#     shell:
#         """
#         salmon quant -i {input.INDEX} -l {params.LIBTYPE} -p {threads} \
#         --validateMappings {params.VALIDATE_MAPPINGS} --gcBias {params.GC_BIAS} \
#         --numGibbsSamples {params.NUM_GIBBS_SAMPLES} -o {output.QUANT} \
#         -1 {input.BAM} 2> {log.log}
#         """

#

# rule ont_sort_index_featureCounts:
#     input:
#         BAM="{OUTDIR}/{SAMPLE}/ont/minimap2/{RECIPE}/sorted_bc.bam.featureCounts.bam",
#     output:
#         BAM="{OUTDIR}/{SAMPLE}/ont/minimap2/{RECIPE}/sorted_bc_gn.bam",
#         BAI="{OUTDIR}/{SAMPLE}/ont/minimap2/{RECIPE}/sorted_bc_gn.bam.bai",
#     run:
#         shell(
#             f"""
#             {EXEC['SAMTOOLS']} sort \
#                 -O BAM \
#                 -o {output.BAM} \
#                 {input.BAM}
#             """
#         )

#         shell(
#             f"""
#             {EXEC['SAMTOOLS']} index {output.BAM}
#             """
#         )


rule ont_add_featureCounts_to_bam:
    input:
        BAM="{OUTDIR}/{SAMPLE}/ont/minimap2/{RECIPE}/sorted_bc_corrected.bam",
        BAI="{OUTDIR}/{SAMPLE}/ont/minimap2/{RECIPE}/sorted_bc_corrected.bam.bai",
        TSV="{OUTDIR}/{SAMPLE}/ont/minimap2/{RECIPE}/sorted_bc_corrected.bam.featureCounts",
    output:
        BAM="{OUTDIR}/{SAMPLE}/ont/minimap2/{RECIPE}/sorted_bc_gn.bam",
    params:
        TAG="GN",
    run:
        shell(
            f"""
            python scripts/py/tsv2tag.py \
                {input.BAM} \
                {input.TSV} \
                {output.BAM} \
                --tag {params.TAG}
            """
        )


# rule ont_index_featureCounts_output:
#     input:
#         BAM="{OUTDIR}/{SAMPLE}/ont/minimap2/{RECIPE}/sorted_bc_gn.bam",
#     output:
#         BAI="{OUTDIR}/{SAMPLE}/ont/minimap2/{RECIPE}/sorted_bc_gn.bam.bai",
#     threads: 1
#     run:
#         shell(
#             f"""
#             {EXEC['SAMTOOLS']} index {input.BAM}
#             """
#         )


# Generate count matrix w/ umi-tools
rule ont_umitools_count:
    input:
        BAM="{OUTDIR}/{SAMPLE}/ont/minimap2/{RECIPE}/sorted_bc_gn.bam",
        BAI="{OUTDIR}/{SAMPLE}/ont/minimap2/{RECIPE}/sorted_bc_gn.bam.bai",
    output:
        COUNTS="{OUTDIR}/{SAMPLE}/ont/minimap2/{RECIPE}/raw/umitools_counts.tsv.gz",
    params:
        CELL_TAG="CB",  # uncorrected = CR
        GENE_TAG="GN",  #GN XS
        UMI_TAG="UR",
    log:
        log="{OUTDIR}/{SAMPLE}/ont/minimap2/{RECIPE}/umitools_count.log",
    threads: 1
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


rule ont_counts_to_sparse:
    input:
        COUNTS="{OUTDIR}/{SAMPLE}/ont/minimap2/{RECIPE}/raw/umitools_counts.tsv.gz",
    output:
        BCS="{OUTDIR}/{SAMPLE}/ont/minimap2/{RECIPE}/raw/barcodes.tsv.gz",
        FEATS="{OUTDIR}/{SAMPLE}/ont/minimap2/{RECIPE}/raw/features.tsv.gz",
        COUNTS="{OUTDIR}/{SAMPLE}/ont/minimap2/{RECIPE}/raw/matrix.mtx.gz",
    params:
        OUTDIR=config["OUTDIR"],
    threads: 1
    run:
        output_dir = output.COUNTS.replace("/matrix.mtx.gz", "")
        shell(
            f"""
            mkdir -p {output_dir}
            python scripts/py/long2mtx.py {input.COUNTS} {output_dir}
            """
        )
