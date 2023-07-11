#############################################
## STAR Alignment
#############################################
# Make output directory, align fastqs, and generate raw/filtered feature/cell-barcode matrices
#   Info for STARsolo command line paramaters: https://github.com/alexdobin/STAR/blob/master/docs/STARsolo.md

#TODO: add multiple recipe compatibility (iterate through the list of space-delimited chemistries listed in sample sheet)
rule STARsolo_align:
    input:
        R1_FQ = '{OUTDIR}/{sample}/tmp/{sample}_R1_final.fq.gz',
        # R1_FQ_HardTrim = '{OUTDIR}/{sample}/tmp/{sample}_R1_HardTrim.fq.gz',
        # R1_FQ_InternalTrim = '{OUTDIR}/{sample}/tmp/{sample}_R1_InternalTrim.fq.gz',
        R2_FQ = '{OUTDIR}/{sample}/tmp/{sample}_R2_final.fq.gz',
        R1_FQ_FILTERED = '{OUTDIR}/{sample}/tmp/{sample}_R1_final_filtered.fq.gz',
        R2_FQ_FILTERED = '{OUTDIR}/{sample}/tmp/{sample}_R2_final_filtered.fq.gz',
        BB_WHITELIST = "{OUTDIR}/{sample}/bb/whitelist.txt",
        BB_1 = "{OUTDIR}/{sample}/bb/whitelist_1.txt",
        BB_2 = "{OUTDIR}/{sample}/bb/whitelist_2.txt",
        BB_ADAPTER = "{OUTDIR}/{sample}/bb/whitelist_adapter.txt"
    output:
        SORTEDBAM = '{OUTDIR}/{sample}/STARsolo/Aligned.sortedByCoord.out.bam', #TODO: add temp()
        UNMAPPED1 = '{OUTDIR}/{sample}/STARsolo/Unmapped.out.mate1',
        UNMAPPED2 = '{OUTDIR}/{sample}/STARsolo/Unmapped.out.mate2',
        VEL_MAT = "{OUTDIR}/{sample}/STARsolo/Solo.out/Velocyto/raw/spliced.mtx",
        VEL_BC = "{OUTDIR}/{sample}/STARsolo/Solo.out/Velocyto/raw/barcodes.tsv",
        VEL_FEAT = "{OUTDIR}/{sample}/STARsolo/Solo.out/Velocyto/raw/features.tsv", 
        GENE_MAT = "{OUTDIR}/{sample}/STARsolo/Solo.out/Gene/raw/matrix.mtx", 
        GENE_BC = "{OUTDIR}/{sample}/STARsolo/Solo.out/Gene/raw/barcodes.tsv", 
        GENE_FEAT = "{OUTDIR}/{sample}/STARsolo/Solo.out/Gene/raw/features.tsv",
        GENEFULL_MAT = "{OUTDIR}/{sample}/STARsolo/Solo.out/GeneFull/raw/matrix.mtx",
        GENEFULL_MAT_EM = "{OUTDIR}/{sample}/STARsolo/Solo.out/GeneFull/raw/UniqueAndMult-EM.mtx",
        GENEFULL_BC = "{OUTDIR}/{sample}/STARsolo/Solo.out/GeneFull/raw/barcodes.tsv",
        GENEFULL_FEAT = "{OUTDIR}/{sample}/STARsolo/Solo.out/GeneFull/raw/features.tsv"
    params:
        OUTDIR = config['OUTDIR'],
        MEMLIMIT = config['MEMLIMIT']
    threads:
        config['CORES']
    priority:
        42
    run:
        tmp_recipe = RECIPE_DICT[wildcards.sample]
        STAR_REF = REF_DICT[wildcards.sample] # Use rRNA reference
        nBB = sum(1 for line in open(input.BB_WHITELIST)) # get number of bead barcodes for filtered count matrix, `--soloCellFilter`

        #TODO: add try catches
        soloType = RECIPE_SHEET["STAR.soloType"][tmp_recipe]
        soloUMI = RECIPE_SHEET["STAR.soloUMI"][tmp_recipe]
        soloCB = RECIPE_SHEET["STAR.soloCB"][tmp_recipe]
        soloCBmatchWLtype = RECIPE_SHEET["STAR.soloCBmatchWLtype"][tmp_recipe]
        soloAdapter = RECIPE_SHEET["STAR.soloAdapter"][tmp_recipe]
        extraSTAR = RECIPE_SHEET["STAR.extra"][tmp_recipe]

        #param handling for different alignment strategies
        if "stomics" in tmp_recipe:
            whitelist = input.BB_WHITELIST
        elif "noTrim" in tmp_recipe:
            # ["seeker_v3.1_noTrimMatchLinker","seeker_v3.1_noTrim_total"]:
            whitelist = f"{input.BB_1} {input.BB_2}"
            # R1 = input.R1_FQ
        elif "internalTrim" in tmp_recipe:
            # ["seeker_v3.1_internalTrim_total"]:
            whitelist = input.BB_WHITELIST
            # R1 = input.R1_FQ_InternalTrim
        elif "adapterInsert" in tmp_recipe:
            whitelist = input.BB_ADAPTER
            # R1 = input.R1_FQ
        else:
            whitelist = input.BB_WHITELIST
            # R1 = input.R1_FQ_HardTrim

        # Select R2 based on alignment recipe
        if "rRNA" in tmp_recipe: # Use trimmed & rRNA-filtered .fq's
            R1 = input.R1_FQ_FILTERED
            R2 = input.R2_FQ_FILTERED
        else: # just trimmed .fq's
            R1 = input.R1_FQ
            R2 = input.R2_FQ

        # Run STARsolo
        shell(
            f"""
            mkdir -p {params.OUTDIR}/{wildcards.sample}/STARsolo

            {STAR_EXEC} \
            --runThreadN {threads} \
            --outFileNamePrefix {params.OUTDIR}/{wildcards.sample}/STARsolo/ \
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
        VEL_MAT = "{OUTDIR}/{sample}/STARsolo/Solo.out/Velocyto/raw/spliced.mtx",
        VEL_BC = "{OUTDIR}/{sample}/STARsolo/Solo.out/Velocyto/raw/barcodes.tsv",
        VEL_FEAT = "{OUTDIR}/{sample}/STARsolo/Solo.out/Velocyto/raw/features.tsv", 
        GENE_MAT = "{OUTDIR}/{sample}/STARsolo/Solo.out/Gene/raw/matrix.mtx", 
        GENE_BC = "{OUTDIR}/{sample}/STARsolo/Solo.out/Gene/raw/barcodes.tsv", 
        GENE_FEAT = "{OUTDIR}/{sample}/STARsolo/Solo.out/Gene/raw/features.tsv",
        GENEFULL_MAT = "{OUTDIR}/{sample}/STARsolo/Solo.out/GeneFull/raw/matrix.mtx",
        GENEFULL_MAT_EM = "{OUTDIR}/{sample}/STARsolo/Solo.out/GeneFull/raw/UniqueAndMult-EM.mtx",
        GENEFULL_BC = "{OUTDIR}/{sample}/STARsolo/Solo.out/GeneFull/raw/barcodes.tsv",
        GENEFULL_FEAT = "{OUTDIR}/{sample}/STARsolo/Solo.out/GeneFull/raw/features.tsv"
    output:
        VEL_MAT = "{OUTDIR}/{sample}/STARsolo/Solo.out/Velocyto/raw/spliced.mtx.gz",
        VEL_BC = "{OUTDIR}/{sample}/STARsolo/Solo.out/Velocyto/raw/barcodes.tsv.gz",
        VEL_FEAT = "{OUTDIR}/{sample}/STARsolo/Solo.out/Velocyto/raw/features.tsv.gz", 
        GENE_MAT = "{OUTDIR}/{sample}/STARsolo/Solo.out/Gene/raw/matrix.mtx.gz", 
        GENE_BC = "{OUTDIR}/{sample}/STARsolo/Solo.out/Gene/raw/barcodes.tsv.gz", 
        GENE_FEAT = "{OUTDIR}/{sample}/STARsolo/Solo.out/Gene/raw/features.tsv.gz",
        GENEFULL_MAT = "{OUTDIR}/{sample}/STARsolo/Solo.out/GeneFull/raw/matrix.mtx.gz",
        GENEFULL_MAT_EM = "{OUTDIR}/{sample}/STARsolo/Solo.out/GeneFull/raw/UniqueAndMult-EM.mtx.gz",
        GENEFULL_BC = "{OUTDIR}/{sample}/STARsolo/Solo.out/GeneFull/raw/barcodes.tsv.gz",
        GENEFULL_FEAT = "{OUTDIR}/{sample}/STARsolo/Solo.out/GeneFull/raw/features.tsv.gz"
    params:
        VELDIR = directory("{OUTDIR}/{sample}/STARsolo/Solo.out/Velocyto"),
        GENEDIR = directory("{OUTDIR}/{sample}/STARsolo/Solo.out/Gene"),
        GENEFULLDIR = directory("{OUTDIR}/{sample}/STARsolo/Solo.out/GeneFull")
    threads:
        config["CORES"]
    run:
        tmp_recipe = RECIPE_DICT[wildcards.sample]
        if "noTrim" in tmp_recipe and "seeker" in tmp_recipe:
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
            pigz -p{threads} -f {OUTDIR}/{wildcards.sample}/STARsolo/*/*/*/*.tsv {OUTDIR}/{wildcards.sample}/STARsolo/*/*/*/*.mtx
            """
        )


# Index .bam output by STAR
rule indexSortedBAM:
    input:
        SORTEDBAM = '{OUTDIR}/{sample}/STARsolo/Aligned.sortedByCoord.out.bam'
    output:
        BAI = '{OUTDIR}/{sample}/STARsolo/Aligned.sortedByCoord.out.bam.bai'
    threads:
        config['CORES']
        # 1
    run:
        shell(
            f"""
            {SAMTOOLS_EXEC} index -@ {threads} {input.SORTEDBAM}
            """
        )