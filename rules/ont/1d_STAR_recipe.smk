
rule ont_STARsolo_align:
    input:
        # R1_FQ = "{OUTDIR}/{SAMPLE}/ont/adapter_scan_readids/full_len_R1.fq.gz",
        # R2_FQ = "{OUTDIR}/{SAMPLE}/ont/adapter_scan_readids/full_len_R2.fq.gz",
        R1_FQ = '{OUTDIR}/{SAMPLE}/ont/tmp/cut_R1.fq.gz',
        R2_FQ = '{OUTDIR}/{SAMPLE}/ont/tmp/cut_R2.fq.gz',
        BB_WHITELIST = "{OUTDIR}/{SAMPLE}/bb/whitelist.txt",
        BB_1 = "{OUTDIR}/{SAMPLE}/bb/whitelist_1.txt",
        BB_2 = "{OUTDIR}/{SAMPLE}/bb/whitelist_2.txt",
        BB_ADAPTER = "{OUTDIR}/{SAMPLE}/bb/whitelist_adapter.txt"
    output:
        BAM = "{OUTDIR}/{SAMPLE}/ont/STARsolo/{RECIPE}/Aligned.sortedByCoord.out.bam", #TODO: add temp()
        UNMAPPED1 = "{OUTDIR}/{SAMPLE}/ont/STARsolo/{RECIPE}/Unmapped.out.mate1",
        UNMAPPED2 = "{OUTDIR}/{SAMPLE}/ont/STARsolo/{RECIPE}/Unmapped.out.mate2",
        VEL_MAT = "{OUTDIR}/{SAMPLE}/ont/STARsolo/{RECIPE}/Solo.out/Velocyto/raw/spliced.mtx",
        VEL_BC = "{OUTDIR}/{SAMPLE}/ont/STARsolo/{RECIPE}/Solo.out/Velocyto/raw/barcodes.tsv",
        VEL_FEAT = "{OUTDIR}/{SAMPLE}/ont/STARsolo/{RECIPE}/Solo.out/Velocyto/raw/features.tsv", 
        GENEFULL_MAT = "{OUTDIR}/{SAMPLE}/ont/STARsolo/{RECIPE}/Solo.out/GeneFull/raw/matrix.mtx",
        GENEFULL_MAT_EM = "{OUTDIR}/{SAMPLE}/ont/STARsolo/{RECIPE}/Solo.out/GeneFull/raw/UniqueAndMult-EM.mtx",
        GENEFULL_BC = "{OUTDIR}/{SAMPLE}/ont/STARsolo/{RECIPE}/Solo.out/GeneFull/raw/barcodes.tsv",
        GENEFULL_FEAT = "{OUTDIR}/{SAMPLE}/ont/STARsolo/{RECIPE}/Solo.out/GeneFull/raw/features.tsv"
    params:
        MEMLIMIT = config["MEMLIMIT"]
    threads:
        config["CORES"]
    resources:
        mem_mb = config["MEMLIMIT_MB"]
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
        elif "noTrim" in recipe:
            whitelist = f"{input.BB_1} {input.BB_2}"
        elif "internalTrim" in recipe:
            whitelist = input.BB_WHITELIST
        elif "adapterInsert" in recipe:
            whitelist = input.BB_ADAPTER
        else:
            whitelist = input.BB_WHITELIST

        # Select input reads based on alignment recipe
        R1 = input.R1_FQ
        R2 = input.R2_FQ

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
#

# compress outputs from STAR (count matrices, cell barcodes, and gene lists)
rule ont_compress_STAR_outs:
    input:
        VEL_MAT = "{OUTDIR}/{SAMPLE}/ont/STARsolo/{RECIPE}/Solo.out/Velocyto/raw/spliced.mtx",
        VEL_BC = "{OUTDIR}/{SAMPLE}/ont/STARsolo/{RECIPE}/Solo.out/Velocyto/raw/barcodes.tsv",
        VEL_FEAT = "{OUTDIR}/{SAMPLE}/ont/STARsolo/{RECIPE}/Solo.out/Velocyto/raw/features.tsv", 
        GENEFULL_MAT = "{OUTDIR}/{SAMPLE}/ont/STARsolo/{RECIPE}/Solo.out/GeneFull/raw/matrix.mtx",
        GENEFULL_MAT_EM = "{OUTDIR}/{SAMPLE}/ont/STARsolo/{RECIPE}/Solo.out/GeneFull/raw/UniqueAndMult-EM.mtx",
        GENEFULL_BC = "{OUTDIR}/{SAMPLE}/ont/STARsolo/{RECIPE}/Solo.out/GeneFull/raw/barcodes.tsv",
        GENEFULL_FEAT = "{OUTDIR}/{SAMPLE}/ont/STARsolo/{RECIPE}/Solo.out/GeneFull/raw/features.tsv"
    output:
        VEL_MAT = "{OUTDIR}/{SAMPLE}/ont/STARsolo/{RECIPE}/Solo.out/Velocyto/raw/spliced.mtx.gz",
        VEL_BC = "{OUTDIR}/{SAMPLE}/ont/STARsolo/{RECIPE}/Solo.out/Velocyto/raw/barcodes.tsv.gz",
        VEL_FEAT = "{OUTDIR}/{SAMPLE}/ont/STARsolo/{RECIPE}/Solo.out/Velocyto/raw/features.tsv.gz",
        GENEFULL_MAT = "{OUTDIR}/{SAMPLE}/ont/STARsolo/{RECIPE}/Solo.out/GeneFull/raw/matrix.mtx.gz",
        GENEFULL_MAT_EM = "{OUTDIR}/{SAMPLE}/ont/STARsolo/{RECIPE}/Solo.out/GeneFull/raw/UniqueAndMult-EM.mtx.gz",
        GENEFULL_BC = "{OUTDIR}/{SAMPLE}/ont/STARsolo/{RECIPE}/Solo.out/GeneFull/raw/barcodes.tsv.gz",
        GENEFULL_FEAT = "{OUTDIR}/{SAMPLE}/ont/STARsolo/{RECIPE}/Solo.out/GeneFull/raw/features.tsv.gz"
    params:
        VELDIR = directory("{OUTDIR}/{SAMPLE}/ont/STARsolo/{RECIPE}/Solo.out/Velocyto"),
        GENEFULLDIR = directory("{OUTDIR}/{SAMPLE}/ont/STARsolo/{RECIPE}/Solo.out/GeneFull")
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

                cat {params.GENEFULLDIR}/raw/barcodes.tsv | sed 's/_//' > {params.GENEFULLDIR}/raw/barcodes_noUnderscore.tsv
                cat {params.GENEFULLDIR}/filtered/barcodes.tsv | sed 's/_//' > {params.GENEFULLDIR}/filtered/barcodes_noUnderscore.tsv
                """
            )

        shell(
            f"""
            {EXEC['PIGZ']} -p{threads} -f \
                {OUTDIR}/{wildcards.SAMPLE}/ont/STARsolo/{recipe}/*/*/*/*.tsv \
                {OUTDIR}/{wildcards.SAMPLE}/ont/STARsolo/{recipe}/*/*/*/*.mtx
            """
        )
#

# Index .bam output by STAR
rule ont_indexSortedBAM:
    input:
        BAM = "{OUTDIR}/{SAMPLE}/ont/STARsolo/{RECIPE}/Aligned.sortedByCoord.out.bam"
    output:
        BAI = "{OUTDIR}/{SAMPLE}/ont/STARsolo/{RECIPE}/Aligned.sortedByCoord.out.bam.bai"
    threads:
        config["CORES"]
        # 1
    run:
        shell(
            f"""
            {EXEC['SAMTOOLS']} index -@ {threads} {input.BAM}
            """
        )
#