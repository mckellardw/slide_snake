#############################################
## STAR Alignment
#############################################
# Make output directory, align fastqs, and generate raw/filtered feature/cell-barcode matrices
#   Info for STARsolo command line paramaters: https://github.com/alexdobin/STAR/blob/master/docs/STARsolo.md

#TODO: add multiple chemistry compatibility (iterate through the list of space-delimited chemistries listed in sample sheet)
rule STARsolo_align:
    input:
        R1_FQ_HardTrim = '{OUTDIR}/{sample}/tmp/{sample}_R1_final.fq.gz',
        R1_FQ = '{OUTDIR}/{sample}/tmp/{sample}_R1_adapterTrim.fq.gz',
        R2_FQ = '{OUTDIR}/{sample}/tmp/{sample}_R2_final.fq.gz',
        FILTERED_R1_FQ = '{OUTDIR}/{sample}/tmp/{sample}_R1_final_filtered.fq.gz',
        FILTERED_R2_FQ = '{OUTDIR}/{sample}/tmp/{sample}_R2_final_filtered.fq.gz',
        BB_WHITELIST = "{OUTDIR}/{sample}/bb/whitelist.txt",
        BB_1 = "{OUTDIR}/{sample}/bb/whitelist_1.txt",
        BB_2 = "{OUTDIR}/{sample}/bb/whitelist_2.txt"
    output:
        SORTEDBAM = '{OUTDIR}/{sample}/STARsolo/Aligned.sortedByCoord.out.bam', #TODO: add temp()
        UNMAPPED1 = '{OUTDIR}/{sample}/STARsolo/Unmapped.out.mate1',
        UNMAPPED2 = '{OUTDIR}/{sample}/STARsolo/Unmapped.out.mate2',
        GENE = directory('{OUTDIR}/{sample}/STARsolo/Solo.out/Gene'),
        GENEFULL = directory('{OUTDIR}/{sample}/STARsolo/Solo.out/GeneFull'),
        # SJ = directory('{OUTDIR}/{sample}/STARsolo/Solo.out/SJ'),
        VEL = directory('{OUTDIR}/{sample}/STARsolo/Solo.out/Velocyto'),
        GENEMAT = '{OUTDIR}/{sample}/STARsolo/Solo.out/Gene/raw/matrix.mtx',
        GENEFULLMAT = '{OUTDIR}/{sample}/STARsolo/Solo.out/GeneFull/raw/matrix.mtx',
        # SJMAT = '{OUTDIR}/{sample}/STARsolo/Solo.out/SJ/raw/matrix.mtx',
        VELMAT = '{OUTDIR}/{sample}/STARsolo/Solo.out/Velocyto/raw/spliced.mtx'
    params:
        OUTDIR = config['OUTDIR'],
        STAR_EXEC = config['STAR_EXEC'],
        MEMLIMIT = config['MEMLIMIT']
    threads:
        config['CORES']
    priority:
        42
    run:
        # print(CHEM_DICT)
        tmp_chemistry = CHEM_DICT[wildcards.sample]
        STAR_REF = REF_DICT[wildcards.sample]
        # BB_WHITELIST = f"{input.BB_1} {input.BB_2}"
        nBB = sum(1 for line in open(input.BB_WHITELIST)) # get number of bead barcodes for filtered count matrix, `--soloCellFilter`

        #TODO: add try catches
        soloType = CHEMISTRY_SHEET["STAR.soloType"][tmp_chemistry]
        soloUMI = CHEMISTRY_SHEET["STAR.soloUMI"][tmp_chemistry]
        soloCB = CHEMISTRY_SHEET["STAR.soloCB"][tmp_chemistry]
        soloCBmatchWLtype = CHEMISTRY_SHEET["STAR.soloCBmatchWLtype"][tmp_chemistry]
        soloAdapter = CHEMISTRY_SHEET["STAR.soloAdapter"][tmp_chemistry]
        extraSTAR = CHEMISTRY_SHEET["STAR.extra"][tmp_chemistry]

        #param handling for different alignment strategies
        if "noTrim" in tmp_chemistry:
            # ["seeker_v3.1_noTrimMatchLinker","seeker_v3.1_noTrim_total"]:
            whitelist = f"{input.BB_1} {input.BB_2}"
            R1 = input.R1_FQ
        else:
            whitelist = input.BB_WHITELIST
            R1 = input.R1_FQ_HardTrim


        if "total" in tmp_chemistry:
            R2 = input.FILTERED_R2_FQ
        else:
            R2 = input.R2_FQ

        shell(
            f"""
            mkdir -p {params.OUTDIR}/{wildcards.sample}/STARsolo

            {params.STAR_EXEC} \
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

            # --outFilterMismatchNoverLmax 0.05 \
            # --outFilterMatchNmin 12 \
            # --outFilterScoreMinOverLread 0 \
            # --outFilterMatchNminOverLread 0 \
        # --soloCellFilter CellRanger2.2 {nBB} 0.99 10 45000 90000 1 0.01 20000 0.01 10000 \

# compress outputs from STAR (count matrices, cell barcodes, and gene lists)
rule compress_STAR_outs:
    input:
        VELMAT = "{OUTDIR}/{sample}/STARsolo/Solo.out/Velocyto/raw/spliced.mtx",
        GENEMAT = "{OUTDIR}/{sample}/STARsolo/Solo.out/Gene/raw/matrix.mtx",
        GENEFULLMAT = "{OUTDIR}/{sample}/STARsolo/Solo.out/GeneFull/raw/matrix.mtx"
    output:
        VELMAT = "{OUTDIR}/{sample}/STARsolo/Solo.out/Velocyto/raw/spliced.mtx.gz",
        GENEMAT = "{OUTDIR}/{sample}/STARsolo/Solo.out/Gene/raw/matrix.mtx.gz",
        GENEFULLMAT = "{OUTDIR}/{sample}/STARsolo/Solo.out/GeneFull/raw/matrix.mtx.gz"
    params:
        VELDIR = directory("{OUTDIR}/{sample}/STARsolo/Solo.out/Velocyto"),
        GENEDIR = directory("{OUTDIR}/{sample}/STARsolo/Solo.out/Gene"),
        GENEFULLDIR = directory("{OUTDIR}/{sample}/STARsolo/Solo.out/GeneFull")
    threads:
        config["CORES_LO"]
    run:
        tmp_chemistry = CHEM_DICT[wildcards.sample]
        if tmp_chemistry in ["seeker_v3.1_noTrimMatchLinker","seeker_v3.1_noTrim_total"]:
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
            pigz -p{threads} {params.VELDIR}/*/*.tsv {params.VELDIR}/*/*.mtx  {params.GENEDIR}/*/*.tsv {params.GENEDIR}/*/*.mtx {params.GENEFULLDIR}/*/*.tsv {params.GENEFULLDIR}/*/*.mtx
            """
        )


# filter the "raw" count matrix by the whitelist
#TODO
# rule STAR_mtx_whitelist_filter:
#     input:
#         # VELMAT = "{OUTDIR}/{sample}/STARsolo/Solo.out/Velocyto/raw/spliced.mtx.gz",
#         GENEMAT = "{OUTDIR}/{sample}/STARsolo/Solo.out/Gene/raw/matrix.mtx.gz",
#         GENEFULLMAT = "{OUTDIR}/{sample}/STARsolo/Solo.out/GeneFull/raw/matrix.mtx.gz"
#     output:
#         # VELMAT = "{OUTDIR}/{sample}/STARsolo/Solo.out/Velocyto/raw_feature_bc_matrix_h5.h5",
#         GENEMAT = "{OUTDIR}/{sample}/STARsolo/Solo.out/Gene/raw_feature_bc_matrix_h5.h5",
#         GENEFULLMAT = "{OUTDIR}/{sample}/STARsolo/Solo.out/GeneFull/raw_feature_bc_matrix_h5.h5"
#     params:
#         # VELDIR = directory("{OUTDIR}/{sample}/STARsolo/Solo.out/Velocyto"),
#         GENEDIR = directory("{OUTDIR}/{sample}/STARsolo/Solo.out/Gene"),
#         GENEFULLDIR = directory("{OUTDIR}/{sample}/STARsolo/Solo.out/GeneFull")
#     # conda:
#     #     "STARsolo"
#     threads:
#         1
#     run:
        #"""mkdir {OUTDIR}/{sample}/STARsolo/Solo.out/Gene/raw_whitelist"""

        #load raw matrix

        #load whitelist

        #filter matrix

        # save mtx 


#TODO add executable
rule indexSortedBAM:
    input:
        SORTEDBAM = '{OUTDIR}/{sample}/STARsolo/Aligned.sortedByCoord.out.bam'
    output:
        BAI = '{OUTDIR}/{sample}/STARsolo/Aligned.sortedByCoord.out.bam.bai'
    threads:
        config['CORES']
    conda:
        "STARsolo"
    shell:
        """
        samtools index -@ {threads} {input.SORTEDBAM}
        """