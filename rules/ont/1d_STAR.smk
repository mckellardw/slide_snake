# STAR rules for ONT data


rule ont_clipBeforeSTAR:
    input:
        FQ="{OUTDIR}/{SAMPLE}/tmp/ont/adapter_scan_readids/full_len_{READ}.fq.gz",
    output:
        FQ="{OUTDIR}/{SAMPLE}/tmp/ont/adapter_scan_readids/full_len_clipped_{READ}.fq.gz",
    params:
        MAX_LENGTH=649,  # 650 is the default length setting for STARlong
    threads: config["CORES"]
    run:
        shell(
            f"""
            zcat {input.FQ} | \
            awk \
                -v maxLength={params.MAX_LENGTH} \
                -f scripts/awk/fq_clipToNBases.awk \
            > {output.FQ.strip('.gz')}

            pigz -p{threads} {output.FQ.strip('.gz')}
            """
        )


# TODO update to match other STARsolo rule(s)
rule ont_STARsolo_align:
    input:
        # R1_FQ = "{OUTDIR}/{SAMPLE}/tmp/ont/adapter_scan_readids/full_len_R1.fq.gz",
        # R2_FQ = "{OUTDIR}/{SAMPLE}/tmp/ont/adapter_scan_readids/full_len_R2.fq.gz",
        R1_FQ="{OUTDIR}/{SAMPLE}/tmp/ont/adapter_scan_readids/full_len_clipped_R1.fq.gz",
        R2_FQ="{OUTDIR}/{SAMPLE}/tmp/ont/adapter_scan_readids/full_len_clipped_R2.fq.gz",
        # R1_FQ = "{OUTDIR}/{SAMPLE}/tmp/ont/cut_R1.fq.gz",
        # R2_FQ = "{OUTDIR}/{SAMPLE}/tmp/ont/cut_R2.fq.gz",
        BC_WHITELIST="{OUTDIR}/{SAMPLE}/bc/whitelist.txt",
        BC_1="{OUTDIR}/{SAMPLE}/bc/whitelist_1.txt",
        BC_2="{OUTDIR}/{SAMPLE}/bc/whitelist_2.txt",
        BC_ADAPTER="{OUTDIR}/{SAMPLE}/bc/whitelist_adapter.txt",
    output:
        BAM="{OUTDIR}/{SAMPLE}/STARsolo/ont/{RECIPE}/Aligned.sortedByCoord.out.bam",  #TODO: add temp()
        UNMAPPED1="{OUTDIR}/{SAMPLE}/STARsolo/ont/{RECIPE}/Unmapped.out.mate1",
        UNMAPPED2="{OUTDIR}/{SAMPLE}/STARsolo/ont/{RECIPE}/Unmapped.out.mate2",
        VEL_MAT="{OUTDIR}/{SAMPLE}/STARsolo/ont/{RECIPE}/Solo.out/Velocyto/raw/spliced.mtx",
        VEL_BC="{OUTDIR}/{SAMPLE}/STARsolo/ont/{RECIPE}/Solo.out/Velocyto/raw/barcodes.tsv",
        VEL_FEAT="{OUTDIR}/{SAMPLE}/STARsolo/ont/{RECIPE}/Solo.out/Velocyto/raw/features.tsv",
        GENEFULL_MAT="{OUTDIR}/{SAMPLE}/STARsolo/ont/{RECIPE}/Solo.out/GeneFull/raw/matrix.mtx",
        GENEFULL_MAT_EM="{OUTDIR}/{SAMPLE}/STARsolo/ont/{RECIPE}/Solo.out/GeneFull/raw/UniqueAndMult-EM.mtx",
        GENEFULL_BC="{OUTDIR}/{SAMPLE}/STARsolo/ont/{RECIPE}/Solo.out/GeneFull/raw/barcodes.tsv",
        GENEFULL_FEAT="{OUTDIR}/{SAMPLE}/STARsolo/ont/{RECIPE}/Solo.out/GeneFull/raw/features.tsv",
    params:
        MEMLIMIT=config["MEMLIMIT"],
        WHITELIST=lambda w: get_whitelist(w),
    threads: config["CORES"]
    resources:
        mem_mb=config["MEMLIMIT_MB"],
    priority: 42
    run:
        # recipe = RECIPE_ONT_DICT[wildcards.SAMPLE]
        recipe = wildcards.RECIPE
        STAR_REF = REF_DICT[wildcards.SAMPLE]
        # nBB = sum(
        #     1 for line in open(input.BC_WHITELIST)
        # )  # get number of bead barcodes for filtered count matrix, `--soloCellFilter`

        # TODO: add try catches
        soloType = RECIPE_SHEET["STAR.soloType"][recipe]
        soloUMI = RECIPE_SHEET["STAR.soloUMI"][recipe]
        soloCB = RECIPE_SHEET["STAR.soloCB"][recipe]
        soloCBmatchWLtype = RECIPE_SHEET["STAR.soloCBmatchWLtype"][recipe]
        soloAdapter = RECIPE_SHEET["STAR.soloAdapter"][recipe]
        extraSTAR = RECIPE_SHEET["STAR.extra"][recipe]
        # extraSTAR = "--outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0 --outFilterMismatchNmax 100 --seedSearchLmax 30 --seedPerReadNmax 100000 --seedPerWindowNmax 100 --alignTranscriptsPerReadNmax 100000 --alignTranscriptsPerWindowNmax 10000"

        # Select input reads based on alignment recipe
        R1 = input.R1_FQ
        R2 = input.R2_FQ

        # Run STARsolo
        # TODO?: --twopassMode
        # WASP?
        # STAR_EXEC = "~/mambaforge/envs/slide_snake/bin/STARlong-avx2"
        shell(
            f"""
            mkdir -p $(dirname {output.BAM})

            STARlong-avx2 \
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
                --soloCBwhitelist {params.WHITELIST} \
                --soloCBmatchWLtype {soloCBmatchWLtype} \
                --soloCellFilter TopCells $(wc -l {params.WHITELIST}) \
                --soloUMIfiltering MultiGeneUMI CR \
                --soloUMIdedup 1MM_CR \
                --soloBarcodeReadLength 0 \
                --soloFeatures Gene GeneFull Velocyto \
                --soloMultiMappers EM
            """
        )


# STAR/source/STARlong --genomeDir /path/to/genomeDir --genomeLoad LoadAndKeep --readFilesIn /path/to/reads.fastq
# --outFileNamePrefix /path/to/output --outSAMtype BAM SortedByCoordinate
# --outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0 --outFilterMismatchNmax 100 --seedSearchLmax 30
# --seedPerReadNmax 100000 --seedPerWindowNmax 100 --alignTranscriptsPerReadNmax 100000
# --alignTranscriptsPerWindowNmax 10000


# compress outputs from STAR (count matrices, cell barcodes, and gene lists)
rule ont_compress_STAR_outs:
    input:
        VEL_MAT="{OUTDIR}/{SAMPLE}/STARsolo/ont/{RECIPE}/Solo.out/Velocyto/raw/spliced.mtx",
        VEL_BC="{OUTDIR}/{SAMPLE}/STARsolo/ont/{RECIPE}/Solo.out/Velocyto/raw/barcodes.tsv",
        VEL_FEAT="{OUTDIR}/{SAMPLE}/STARsolo/ont/{RECIPE}/Solo.out/Velocyto/raw/features.tsv",
        GENEFULL_MAT="{OUTDIR}/{SAMPLE}/STARsolo/ont/{RECIPE}/Solo.out/GeneFull/raw/matrix.mtx",
        GENEFULL_MAT_EM="{OUTDIR}/{SAMPLE}/STARsolo/ont/{RECIPE}/Solo.out/GeneFull/raw/UniqueAndMult-EM.mtx",
        GENEFULL_BC="{OUTDIR}/{SAMPLE}/STARsolo/ont/{RECIPE}/Solo.out/GeneFull/raw/barcodes.tsv",
        GENEFULL_FEAT="{OUTDIR}/{SAMPLE}/STARsolo/ont/{RECIPE}/Solo.out/GeneFull/raw/features.tsv",
    output:
        VEL_MAT="{OUTDIR}/{SAMPLE}/STARsolo/ont/{RECIPE}/Solo.out/Velocyto/raw/spliced.mtx.gz",
        VEL_BC="{OUTDIR}/{SAMPLE}/STARsolo/ont/{RECIPE}/Solo.out/Velocyto/raw/barcodes.tsv.gz",
        VEL_FEAT="{OUTDIR}/{SAMPLE}/STARsolo/ont/{RECIPE}/Solo.out/Velocyto/raw/features.tsv.gz",
        GENEFULL_MAT="{OUTDIR}/{SAMPLE}/STARsolo/ont/{RECIPE}/Solo.out/GeneFull/raw/matrix.mtx.gz",
        GENEFULL_MAT_EM="{OUTDIR}/{SAMPLE}/STARsolo/ont/{RECIPE}/Solo.out/GeneFull/raw/UniqueAndMult-EM.mtx.gz",
        GENEFULL_BC="{OUTDIR}/{SAMPLE}/STARsolo/ont/{RECIPE}/Solo.out/GeneFull/raw/barcodes.tsv.gz",
        GENEFULL_FEAT="{OUTDIR}/{SAMPLE}/STARsolo/ont/{RECIPE}/Solo.out/GeneFull/raw/features.tsv.gz",
    params:
        VELDIR=directory("{OUTDIR}/{SAMPLE}/STARsolo/ont/{RECIPE}/Solo.out/Velocyto"),
        GENEFULLDIR=directory(
            "{OUTDIR}/{SAMPLE}/STARsolo/ont/{RECIPE}/Solo.out/GeneFull"
        ),
    threads: config["CORES"]
    run:
        # recipe = RECIPE_DICT[wildcards.SAMPLE]
        recipe = wildcards.RECIPE

        if "noTrim" in recipe and "seeker" in recipe:
            # ["seeker_v3.1_noTrimMatchLinker","seeker_v3.1_noTrim_total"]:
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
            pigz -p{threads} -f \
                {OUTDIR}/{wildcards.SAMPLE}/STARsolo/ont/{recipe}/*/*/*/*.tsv \
                {OUTDIR}/{wildcards.SAMPLE}/STARsolo/ont/{recipe}/*/*/*/*.mtx
            """
        )
