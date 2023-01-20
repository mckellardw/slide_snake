#############################################
## STAR Alignment
#############################################
# Make output directory, align fastqs, and generate raw/filtered feature/cell-barcode matrices
#   Info for STARsolo command line paramaters: https://github.com/alexdobin/STAR/blob/master/docs/STARsolo.md

rule STARsolo_align:
    input:
        FINAL_R1_FQ = '{OUTDIR}/{sample}/tmp/{sample}_R1_adapterTrim.fq.gz', #*Note*- change to "_R1_Final..." to include adapter hard trimming
        FINAL_R2_FQ = '{OUTDIR}/{sample}/tmp/{sample}_R2_final.fq.gz',
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

        #TODO: check for empty values
        soloType = CHEMISTRY_SHEET["STAR.soloType"][tmp_chemistry]
        soloUMI = CHEMISTRY_SHEET["STAR.soloUMI"][tmp_chemistry]
        soloCB = CHEMISTRY_SHEET["STAR.soloCB"][tmp_chemistry]
        soloCBmatchWLtype = CHEMISTRY_SHEET["STAR.soloCBmatchWLtype"][tmp_chemistry]
        soloAdapter = CHEMISTRY_SHEET["STAR.soloAdapter"][tmp_chemistry]

        if tmp_chemistry == "seeker_v3.1_noTrimMatchLinker":
            whitelist = f"{input.BB_1} {input.BB_2}"
        else:
            whitelist = input.BB_WHITELIST

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
            --readFilesIn {input.FINAL_R2_FQ} {input.FINAL_R1_FQ} \
            --clipAdapterType CellRanger4 \
            --outReadsUnmapped Fastx \
            --outFilterMultimapNmax 50 \
            --outFilterMismatchNoverLmax 0.05 \
            --outFilterMatchNmin 12 \
            --outFilterScoreMinOverLread 0 \
            --outFilterMatchNminOverLread 0 \
            --soloType {soloType} \
            {soloUMI} {soloCB} {soloAdapter} \
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
        1
        # config["CORES_LO"]
    # conda:
    #     "STARsolo"
    run:
        shell(
            f"""
            gzip -qf {params.VELDIR}/*/*.tsv {params.VELDIR}/*/*.mtx
            gzip -qf {params.GENEDIR}/*/*.tsv {params.GENEDIR}/*/*.mtx
            gzip -qf {params.GENEFULLDIR}/*/*.tsv {params.GENEFULLDIR}/*/*.mtx
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
