#############################################
## STAR Alignment
#############################################
# Make output directory, align fastqs, and generate raw/filtered feature/cell-barcode matrices
#   Info for STARsolo command line paramaters: https://github.com/alexdobin/STAR/blob/master/docs/STARsolo.md
rule STARsolo_align:
    input:
        # R1_FQ_HardTrim = '{OUTDIR}/{SAMPLE}/tmp/{SAMPLE}_R1_HardTrim.fq.gz',
        # R1_FQ_InternalTrim = '{OUTDIR}/{SAMPLE}/tmp/{SAMPLE}_R1_InternalTrim.fq.gz',
        R1_FQ = '{OUTDIR}/{SAMPLE}/tmp/cut_R1.fq.gz',
        R2_FQ = '{OUTDIR}/{SAMPLE}/tmp/cut_R2.fq.gz',
        R1_FQ_TWICE_CUT = '{OUTDIR}/{SAMPLE}/tmp/twiceCut_R1.fq.gz',
        R2_FQ_TWICE_CUT = '{OUTDIR}/{SAMPLE}/tmp/twiceCut_R2.fq.gz',
        R1_FQ_STAR_FILTERED = '{OUTDIR}/{SAMPLE}/rRNA/STARsolo/final_filtered_R1.fq.gz',
        R2_FQ_STAR_FILTERED = '{OUTDIR}/{SAMPLE}/rRNA/STARsolo/final_filtered_R2.fq.gz',
        R1_FQ_BWA_FILTERED  = '{OUTDIR}/{SAMPLE}/rRNA/bwa/final_filtered_R1.fq.gz',
        R2_FQ_BWA_FILTERED  = '{OUTDIR}/{SAMPLE}/rRNA/bwa/final_filtered_R2.fq.gz',
        BB_WHITELIST = "{OUTDIR}/{SAMPLE}/bb/whitelist.txt",
        BB_1 = "{OUTDIR}/{SAMPLE}/bb/whitelist_1.txt",
        BB_2 = "{OUTDIR}/{SAMPLE}/bb/whitelist_2.txt",
        BB_ADAPTER = "{OUTDIR}/{SAMPLE}/bb/whitelist_adapter.txt"
    output:
        BAM = '{OUTDIR}/{SAMPLE}/STARsolo/short_read/{RECIPE}/Aligned.sortedByCoord.out.bam', #TODO: add temp()
        UNMAPPED1 = '{OUTDIR}/{SAMPLE}/STARsolo/short_read/{RECIPE}/Unmapped.out.mate1',
        UNMAPPED2 = '{OUTDIR}/{SAMPLE}/STARsolo/short_read/{RECIPE}/Unmapped.out.mate2',
        VEL_MAT = "{OUTDIR}/{SAMPLE}/STARsolo/short_read/{RECIPE}/Solo.out/Velocyto/raw/spliced.mtx",
        VEL_BC = "{OUTDIR}/{SAMPLE}/STARsolo/short_read/{RECIPE}/Solo.out/Velocyto/raw/barcodes.tsv",
        VEL_FEAT = "{OUTDIR}/{SAMPLE}/STARsolo/short_read/{RECIPE}/Solo.out/Velocyto/raw/features.tsv", 
        GENE_MAT = "{OUTDIR}/{SAMPLE}/STARsolo/short_read/{RECIPE}/Solo.out/Gene/raw/matrix.mtx",
        GENE_MAT_EM = "{OUTDIR}/{SAMPLE}/STARsolo/short_read/{RECIPE}/Solo.out/Gene/raw/UniqueAndMult-EM.mtx",
        GENE_BC = "{OUTDIR}/{SAMPLE}/STARsolo/short_read/{RECIPE}/Solo.out/Gene/raw/barcodes.tsv", 
        GENE_FEAT = "{OUTDIR}/{SAMPLE}/STARsolo/short_read/{RECIPE}/Solo.out/Gene/raw/features.tsv",
        GENEFULL_MAT = "{OUTDIR}/{SAMPLE}/STARsolo/short_read/{RECIPE}/Solo.out/GeneFull/raw/matrix.mtx",
        GENEFULL_MAT_EM = "{OUTDIR}/{SAMPLE}/STARsolo/short_read/{RECIPE}/Solo.out/GeneFull/raw/UniqueAndMult-EM.mtx",
        GENEFULL_BC = "{OUTDIR}/{SAMPLE}/STARsolo/short_read/{RECIPE}/Solo.out/GeneFull/raw/barcodes.tsv",
        GENEFULL_FEAT = "{OUTDIR}/{SAMPLE}/STARsolo/short_read/{RECIPE}/Solo.out/GeneFull/raw/features.tsv"# MATS = [f"{OUTDIR}/{SAMPLE}/STARsolo/short_read/{RECIPE}/Solo.out/{SOLO}/raw/{MAT}.mtx" for SAMPLE in SAMPLES for RECIPE in RECIPE_DICT[SAMPLE] for SOLO in ["Gene","GeneFull"] for MAT in ["matrix", "UniqueAndMult-EM"]],
        # VELO_MATS = [f"{OUTDIR}/{SAMPLE}/STARsolo/short_read/{RECIPE}/Solo.out/Velocyto/raw/{MAT}.mtx" for SAMPLE in SAMPLES for RECIPE in RECIPE_DICT[SAMPLE] for MAT in ["spliced","unspliced","ambiguous"]],
        # BC = [f"{OUTDIR}/{SAMPLE}/STARsolo/short_read/{RECIPE}/Solo.out/{SOLO}/raw/barcodes.tsv" for SAMPLE in SAMPLES for RECIPE in RECIPE_DICT[SAMPLE] for SOLO in ["Velocyto","Gene","GeneFull"]],
        # FEAT = [f"{OUTDIR}/{SAMPLE}/STARsolo/short_read/{RECIPE}/Solo.out/{SOLO}/raw/features.tsv" for SAMPLE in SAMPLES for RECIPE in RECIPE_DICT[SAMPLE] for SOLO in ["Velocyto","Gene","GeneFull"]]
    params:
        MEMLIMIT = config['MEMLIMIT'],
    threads:
        config['CORES'],
    # resources:
    #     mem_mb = config['MEMLIMIT_MB'],
    priority:
        42
    run:
        # recipe = RECIPE_DICT[wildcards.SAMPLE]
        recipe = wildcards.RECIPE
        STAR_REF = REF_DICT[wildcards.SAMPLE] # Use rRNA reference
        nBB = sum(1 for line in open(input.BB_WHITELIST)) # get number of bead barcodes for filtered count matrix, `--soloCellFilter`

        #TODO: add try catches
        soloType = RECIPE_SHEET["STAR.soloType"][recipe]
        soloUMI = RECIPE_SHEET["STAR.soloUMI"][recipe]
        soloCB = RECIPE_SHEET["STAR.soloCB"][recipe]
        soloCBmatchWLtype = RECIPE_SHEET["STAR.soloCBmatchWLtype"][recipe]
        soloAdapter = RECIPE_SHEET["STAR.soloAdapter"][recipe]
        extraSTAR = RECIPE_SHEET["STAR.extra"][recipe]

        #param handling for different SlideSeq R1 strategies
        if "stomics" in recipe:
            whitelist = input.BB_WHITELIST
        elif "matchLinker" in recipe:
            whitelist = f"{input.BB_1} {input.BB_2}"
        elif "internalTrim" in recipe:
            whitelist = input.BB_WHITELIST
        elif "adapterInsert" in recipe:
            whitelist = input.BB_ADAPTER
        else:
            whitelist = input.BB_WHITELIST

        # Select input reads based on alignment recipe
        if "rRNA.STAR" in recipe: # Use trimmed & STAR-rRNA-filtered .fq's
            R1 = input.R1_FQ_STAR_FILTERED
            R2 = input.R2_FQ_STAR_FILTERED
        elif "rRNA.bwa" in recipe: #TODO Use trimmed & bwa-rRNA-filtered .fq's
            R1 = input.R1_FQ_BWA_FILTERED
            R2 = input.R2_FQ_BWA_FILTERED
        elif "rRNA" not in recipe: # just trimmed .fq's
            # R1 = input.R1_FQ
            # R2 = input.R2_FQ
            R1 = input.R1_FQ_TWICE_CUT
            R2 = input.R2_FQ_TWICE_CUT
        else:
            print("I just don't know what to do with myself...")

        # Run STARsolo
        #TODO?: --twopassMode
        #WASP?
        shell(
            f"""
            mkdir -p $(dirname {output.BAM})

            {EXEC['STAR']} \
                --runThreadN {threads} \
                --outFileNamePrefix $(dirname {output.BAM})/ \
                --outSAMtype BAM SortedByCoordinate \
                --outSAMattributes NH HI nM AS CR UR CB UB GX GN sS sQ sM \
                --readFilesCommand zcat \
                --genomeDir {STAR_REF} \
                --limitBAMsortRAM={params.MEMLIMIT} \
                --readFilesIn {R2} {R1} \
                --clipAdapterType CellRanger4 \
                --outReadsUnmapped Fastx \
                --outSAMunmapped Within KeepPairs \
                --soloType {soloType} {soloUMI} {soloCB} {soloAdapter} {extraSTAR} \
                --soloCBwhitelist {whitelist} \
                --soloCBmatchWLtype {soloCBmatchWLtype} \
                --soloCellFilter TopCells {nBB} \
                --soloUMIfiltering MultiGeneUMI CR \
                --soloUMIdedup 1MM_CR \
                --soloBarcodeReadLength 0 \
                --soloFeatures Gene GeneFull Velocyto \
                --soloMultiMappers EM
            """
        )

# compress outputs from STAR (count matrices, cell barcodes, and gene lists)
rule compress_STAR_outs:
    input:
        VEL_MAT = "{OUTDIR}/{SAMPLE}/STARsolo/short_read/{RECIPE}/Solo.out/Velocyto/raw/spliced.mtx",
        VEL_BC = "{OUTDIR}/{SAMPLE}/STARsolo/short_read/{RECIPE}/Solo.out/Velocyto/raw/barcodes.tsv",
        VEL_FEAT = "{OUTDIR}/{SAMPLE}/STARsolo/short_read/{RECIPE}/Solo.out/Velocyto/raw/features.tsv", 
        GENE_MAT = "{OUTDIR}/{SAMPLE}/STARsolo/short_read/{RECIPE}/Solo.out/Gene/raw/matrix.mtx",
        GENE_MAT_EM = "{OUTDIR}/{SAMPLE}/STARsolo/short_read/{RECIPE}/Solo.out/Gene/raw/UniqueAndMult-EM.mtx",
        GENE_BC = "{OUTDIR}/{SAMPLE}/STARsolo/short_read/{RECIPE}/Solo.out/Gene/raw/barcodes.tsv", 
        GENE_FEAT = "{OUTDIR}/{SAMPLE}/STARsolo/short_read/{RECIPE}/Solo.out/Gene/raw/features.tsv",
        GENEFULL_MAT = "{OUTDIR}/{SAMPLE}/STARsolo/short_read/{RECIPE}/Solo.out/GeneFull/raw/matrix.mtx",
        GENEFULL_MAT_EM = "{OUTDIR}/{SAMPLE}/STARsolo/short_read/{RECIPE}/Solo.out/GeneFull/raw/UniqueAndMult-EM.mtx",
        GENEFULL_BC = "{OUTDIR}/{SAMPLE}/STARsolo/short_read/{RECIPE}/Solo.out/GeneFull/raw/barcodes.tsv",
        GENEFULL_FEAT = "{OUTDIR}/{SAMPLE}/STARsolo/short_read/{RECIPE}/Solo.out/GeneFull/raw/features.tsv"
        # MATS = [f"{OUTDIR}/{SAMPLE}/STARsolo/short_read/{RECIPE}/Solo.out/{SOLO}/raw/{MAT}.mtx" for SAMPLE in SAMPLES for RECIPE in RECIPE_DICT[SAMPLE] for SOLO in ["Gene","GeneFull"] for MAT in ["matrix", "UniqueAndMult-EM"]],
        # VELO_MATS = [f"{OUTDIR}/{SAMPLE}/STARsolo/short_read/{RECIPE}/Solo.out/Velocyto/raw/{MAT}.mtx" for SAMPLE in SAMPLES for RECIPE in RECIPE_DICT[SAMPLE] for MAT in ["spliced","unspliced","ambiguous"]],
        # BC = [f"{OUTDIR}/{SAMPLE}/STARsolo/short_read/{RECIPE}/Solo.out/{SOLO}/raw/barcodes.tsv" for SAMPLE in SAMPLES for RECIPE in RECIPE_DICT[SAMPLE] for SOLO in ["Velocyto","Gene","GeneFull"]],
        # FEAT = [f"{OUTDIR}/{SAMPLE}/STARsolo/short_read/{RECIPE}/Solo.out/{SOLO}/raw/features.tsv" for SAMPLE in SAMPLES for RECIPE in RECIPE_DICT[SAMPLE] for SOLO in ["Velocyto","Gene","GeneFull"]]
    output:
        # MATS = [f"{OUTDIR}/{SAMPLE}/STARsolo/short_read/{RECIPE}/Solo.out/{SOLO}/raw/{MAT}.mtx.gz" for SAMPLE in SAMPLES for RECIPE in RECIPE_DICT[SAMPLE] for SOLO in ["Gene","GeneFull"] for MAT in ["matrix", "UniqueAndMult-EM"]],
        # VELO_MATS = [f"{OUTDIR}/{SAMPLE}/STARsolo/short_read/{RECIPE}/Solo.out/Velocyto/raw/{MAT}.mtx.gz" for SAMPLE in SAMPLES for RECIPE in RECIPE_DICT[SAMPLE] for MAT in ["spliced","unspliced","ambiguous"]],
        # BC = [f"{OUTDIR}/{SAMPLE}/STARsolo/short_read/{RECIPE}/Solo.out/{SOLO}/raw/barcodes.tsv.gz" for SAMPLE in SAMPLES for RECIPE in RECIPE_DICT[SAMPLE] for SOLO in ["Velocyto","Gene","GeneFull"]],
        # FEAT = [f"{OUTDIR}/{SAMPLE}/STARsolo/short_read/{RECIPE}/Solo.out/{SOLO}/raw/features.tsv.gz" for SAMPLE in SAMPLES for RECIPE in RECIPE_DICT[SAMPLE] for SOLO in ["Velocyto","Gene","GeneFull"]]
        VEL_MAT = "{OUTDIR}/{SAMPLE}/STARsolo/short_read/{RECIPE}/Solo.out/Velocyto/raw/spliced.mtx.gz",
        VEL_BC = "{OUTDIR}/{SAMPLE}/STARsolo/short_read/{RECIPE}/Solo.out/Velocyto/raw/barcodes.tsv.gz",
        VEL_FEAT = "{OUTDIR}/{SAMPLE}/STARsolo/short_read/{RECIPE}/Solo.out/Velocyto/raw/features.tsv.gz", 
        GENE_MAT = "{OUTDIR}/{SAMPLE}/STARsolo/short_read/{RECIPE}/Solo.out/Gene/raw/matrix.mtx.gz",
        GENE_MAT_EM = "{OUTDIR}/{SAMPLE}/STARsolo/short_read/{RECIPE}/Solo.out/Gene/raw/UniqueAndMult-EM.mtx.gz",
        GENE_BC = "{OUTDIR}/{SAMPLE}/STARsolo/short_read/{RECIPE}/Solo.out/Gene/raw/barcodes.tsv.gz", 
        GENE_FEAT = "{OUTDIR}/{SAMPLE}/STARsolo/short_read/{RECIPE}/Solo.out/Gene/raw/features.tsv.gz",
        GENEFULL_MAT = "{OUTDIR}/{SAMPLE}/STARsolo/short_read/{RECIPE}/Solo.out/GeneFull/raw/matrix.mtx.gz",
        GENEFULL_MAT_EM = "{OUTDIR}/{SAMPLE}/STARsolo/short_read/{RECIPE}/Solo.out/GeneFull/raw/UniqueAndMult-EM.mtx.gz",
        GENEFULL_BC = "{OUTDIR}/{SAMPLE}/STARsolo/short_read/{RECIPE}/Solo.out/GeneFull/raw/barcodes.tsv.gz",
        GENEFULL_FEAT = "{OUTDIR}/{SAMPLE}/STARsolo/short_read/{RECIPE}/Solo.out/GeneFull/raw/features.tsv.gz"
    params:
        VELDIR = directory("{OUTDIR}/{SAMPLE}/STARsolo/short_read/{RECIPE}/Solo.out/Velocyto"),
        GENEDIR = directory("{OUTDIR}/{SAMPLE}/STARsolo/short_read/{RECIPE}/Solo.out/Gene"),
        GENEFULLDIR = directory("{OUTDIR}/{SAMPLE}/STARsolo/short_read/{RECIPE}/Solo.out/GeneFull")
    threads:
        config["CORES"]
    run:
        # recipe = RECIPE_DICT[wildcards.SAMPLE]
        recipe = wildcards.RECIPE

        if "noTrim" in recipe and "seeker" in recipe:
            #["seeker_v3.1_noTrimMatchLinker","seeker_v3.1_noTrim_total"]:
            shell(
                f"""
                cat {params.VELDIR}/raw/barcodes.tsv | sed 's/_//' > {params.VELDIR}/raw/barcodes_noUnderscore.tsv
                cat {params.VELDIR}/filtered/barcodes.tsv | sed 's/_//' > {params.VELDIR}/filtered/barcodes_noUnderscore.tsv

                cat {params.GENEDIR}/raw/barcodes.tsv | sed 's/_//' > {params.GENEDIR}/raw/barcodes_noUnderscore.tsv
                cat {params.GENEDIR}/filtered/barcodes.tsv | sed 's/_//' > {params.GENEDIR}/filtered/barcodes_noUnderscore.tsv

                cat {params.GENEFULLDIR}/raw/barcodes.tsv | sed 's/_//' > {params.GENEFULLDIR}/raw/barcodes_noUnderscore.tsv
                cat {params.GENEFULLDIR}/filtered/barcodes.tsv | sed 's/_//' > {params.GENEFULLDIR}/filtered/barcodes_noUnderscore.tsv
                """
            )

        shell(
            f"""
            {EXEC['PIGZ']} -p{threads} -f \
                {OUTDIR}/{wildcards.SAMPLE}/STARsolo/short_read/{recipe}/*/*/*/*.tsv \
                {OUTDIR}/{wildcards.SAMPLE}/STARsolo/short_read/{recipe}/*/*/*/*.mtx
            """
        )


# Index .bam output by STAR
rule indexSortedBAM:
    input:
        BAM = '{OUTDIR}/{SAMPLE}/STARsolo/short_read/{RECIPE}/Aligned.sortedByCoord.out.bam'
    output:
        BAI = '{OUTDIR}/{SAMPLE}/STARsolo/short_read/{RECIPE}/Aligned.sortedByCoord.out.bam.bai'
    threads:
        config['CORES']
        # 1
    run:
        shell(
            f"""
            {EXEC['SAMTOOLS']} index -@ {threads} {input.BAM}
            """
        )


# rule multiqc_STAR:
#     input:
#         expand(
#             "{OUTDIR}/{SAMPLE}/STARsolo/short_read/{RECIPE}/Aligned.sortedByCoord.out.bam", 
#             SAMPLE=SAMPLES, RECIPE=RECIPES),
#         expand("{OUTDIR}/{SAMPLE}/STARsolo/short_read/{RECIPE}/Unmapped.out.mate1", 
#                SAMPLE=SAMPLES, RECIPE=RECIPES),
#         expand("{OUTDIR}/{SAMPLE}/STARsolo/short_read/{RECIPE}/Unmapped.out.mate2", 
#                SAMPLE=SAMPLES, RECIPE=RECIPES),
#         expand("{OUTDIR}/{SAMPLE}/STARsolo/short_read/{RECIPE}/Solo.out/Velocyto/raw/spliced.mtx", 
#                SAMPLE=SAMPLES, RECIPE=RECIPES),
#         expand("{OUTDIR}/{SAMPLE}/STARsolo/short_read/{RECIPE}/Solo.out/Velocyto/raw/barcodes.tsv", 
#                SAMPLE=SAMPLES, RECIPE=RECIPES),
#         expand("{OUTDIR}/{SAMPLE}/STARsolo/short_read/{RECIPE}/Solo.out/Velocyto/raw/features.tsv", 
#                SAMPLE=SAMPLES, RECIPE=RECIPES),
#         expand("{OUTDIR}/{SAMPLE}/STARsolo/short_read/{RECIPE}/Solo.out/Gene/raw/matrix.mtx", 
#                SAMPLE=SAMPLES, RECIPE=RECIPES),
#         expand("{OUTDIR}/{SAMPLE}/STARsolo/short_read/{RECIPE}/Solo.out/Gene/raw/UniqueAndMult-EM.mtx", 
#                SAMPLE=SAMPLES, RECIPE=RECIPES),
#         expand("{OUTDIR}/{SAMPLE}/STARsolo/short_read/{RECIPE}/Solo.out/Gene/raw/barcodes.tsv", 
#                SAMPLE=SAMPLES, RECIPE=RECIPES),
#         expand("{OUTDIR}/{SAMPLE}/STARsolo/short_read/{RECIPE}/Solo.out/Gene/raw/features.tsv", 
#                SAMPLE=SAMPLES, RECIPE=RECIPES),
#         expand("{OUTDIR}/{SAMPLE}/STARsolo/short_read/{RECIPE}/Solo.out/GeneFull/raw/matrix.mtx", 
#                SAMPLE=SAMPLES, RECIPE=RECIPES),
#         expand("{OUTDIR}/{SAMPLE}/STARsolo/short_read/{RECIPE}/Solo.out/GeneFull/raw/UniqueAndMult-EM.mtx", 
#                SAMPLE=SAMPLES, RECIPE=RECIPES),
#         expand("{OUTDIR}/{SAMPLE}/STARsolo/short_read/{RECIPE}/Solo.out/GeneFull/raw/barcodes.tsv", 
#                SAMPLE=SAMPLES, RECIPE=RECIPES),
#         expand("{OUTDIR}/{SAMPLE}/STARsolo/short_read/{RECIPE}/Solo.out/GeneFull/raw/features.tsv", 
#                SAMPLE=SAMPLES, RECIPE=RECIPES)
#     output:
#         REPORT = "{OUTDIR}/{SAMPLE}/STARsolo/short_read/multiqc_report.html"
#     params:
#         config['MEMLIMIT']
#     threads:
#         config['CORES']
#     resources:
#         mem_mb = config['MEMLIMIT_MB']
#     shell:
#         """
#         multiqc . -o {output.REPORT}
#         """
