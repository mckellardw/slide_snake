#############################################
## STAR Alignment
#############################################
# Make output directory, align fastqs, and generate raw/filtered feature/cell-barcode matrices
#   Info for STARsolo command line paramaters: https://github.com/alexdobin/STAR/blob/master/docs/STARsolo.md

#TODO: add multiple chemistry compatibility (iterate through the list of space-delimited chemistries listed in sample sheet)
rule STARsolo_align:
    input:
        R1_FQ = '{OUTDIR}/{sample}/tmp/{sample}_R1_adapterTrim.fq.gz',
        R1_FQ_HardTrim = '{OUTDIR}/{sample}/tmp/{sample}_R1_finalHardTrim.fq.gz',
        R1_FQ_InternalTrim = '{OUTDIR}/{sample}/tmp/{sample}_R1_finalInternalTrim.fq.gz',
        R2_FQ = '{OUTDIR}/{sample}/tmp/{sample}_R2_final.fq.gz',
        R1_FQ_FILTERED = '{OUTDIR}/{sample}/tmp/{sample}_R1_final_filtered.fq.gz',
        R2_FQ_FILTERED = '{OUTDIR}/{sample}/tmp/{sample}_R2_final_filtered.fq.gz',
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
        # STAR_EXEC = config['STAR_EXEC'],
        MEMLIMIT = config['MEMLIMIT']
    threads:
        config['CORES']
    priority:
        42
    run:
        tmp_chemistry = CHEM_DICT[wildcards.sample]
        STAR_REF = REF_DICT[wildcards.sample] # Use rRNA reference
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
        elif "internalTrim" in tmp_chemistry:
            # ["seeker_v3.1_internalTrim_total"]:
            whitelist = input.BB_WHITELIST
            R1 = input.R1_FQ_InternalTrim
        else:
            whitelist = input.BB_WHITELIST
            R1 = input.R1_FQ_HardTrim

        # Select R2 based on alignment recipe
        if "total" in tmp_chemistry:
            R1 = input.R1_FQ_FILTERED
            R2 = input.R2_FQ_FILTERED
        else:
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
        config["CORES"]
    run:
        tmp_chemistry = CHEM_DICT[wildcards.sample]
        if "noTrim" in tmp_chemistry:
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
            pigz -p{threads} {OUTDIR}/{wildcards.sample}/STARsolo/*/*/*/*.tsv {OUTDIR}/{wildcards.sample}/STARsolo/*/*/*/*.mtx
            """
        )


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
        {SAMTOOLS_EXEC} index -@ {threads} {input.SORTEDBAM}
        """
