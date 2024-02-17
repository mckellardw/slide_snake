# Make output directory, align fastqs, and generate raw/filtered feature/cell-barcode matrices
#   Info for STARsolo command line paramaters: https://github.com/alexdobin/STAR/blob/master/docs/STARsolo.md
#TODO- dedup rRNA .bam files
rule STAR_rRNA_align:
    input:
        R1_FQ = '{OUTDIR}/{SAMPLE}/tmp/twiceCut_R1.fq.gz',
        R2_FQ = '{OUTDIR}/{SAMPLE}/tmp/twiceCut_R2.fq.gz',
        BB_WHITELIST = "{OUTDIR}/{SAMPLE}/bb/whitelist.txt",
        BB_1 = "{OUTDIR}/{SAMPLE}/bb/whitelist_1.txt",
        BB_2 = "{OUTDIR}/{SAMPLE}/bb/whitelist_2.txt",
        BB_ADAPTER = "{OUTDIR}/{SAMPLE}/bb/whitelist_adapter.txt"
    output:
        BAM = '{OUTDIR}/{SAMPLE}/rRNA/STARsolo/Aligned.sortedByCoord.out.bam', #TODO: add temp()
        UNMAPPED1 = '{OUTDIR}/{SAMPLE}/rRNA/STARsolo/Unmapped.out.mate1',
        UNMAPPED2 = '{OUTDIR}/{SAMPLE}/rRNA/STARsolo/Unmapped.out.mate2',
        GENEDIRECTORY = directory('{OUTDIR}/{SAMPLE}/rRNA/STARsolo/Solo.out/GeneFull'),
        GENEMAT = '{OUTDIR}/{SAMPLE}/rRNA/STARsolo/Solo.out/GeneFull/raw/matrix.mtx'
    params:
        MEMLIMIT = config['MEMLIMIT']
    threads:
        config['CORES']
    run: 
        recipe = RECIPE_DICT[wildcards.SAMPLE][0] #TODO fix recipe handling here...
        # recipe = wildcards.RECIPE 
        
        STAR_REF = rRNA_STAR_DICT[wildcards.SAMPLE] # use rRNA ref
        nBB = sum(1 for line in open(input.BB_WHITELIST)) # get number of bead barcodes for filtered count matrix, `--soloCellFilter`

        #TODO: add try catches
        soloType = RECIPE_SHEET["STAR.soloType"][recipe]
        soloUMI = RECIPE_SHEET["STAR.soloUMI"][recipe]
        soloCB = RECIPE_SHEET["STAR.soloCB"][recipe]
        soloCBmatchWLtype = RECIPE_SHEET["STAR.soloCBmatchWLtype"][recipe]
        soloAdapter = RECIPE_SHEET["STAR.soloAdapter"][recipe]
        extraSTAR = RECIPE_SHEET["STAR.extra"][recipe]

        # Params for different barcode handling strategies
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
                --readFilesIn {input.R2_FQ} {input.R1_FQ} \
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
                --soloFeatures GeneFull \
                --soloMultiMappers EM
            """
        )

# compress outputs from STAR (count matrices, cell barcodes, and gene lists)
rule STAR_rRNA_compress_outs:
    input:
        GENEMAT = "{OUTDIR}/{SAMPLE}/rRNA/STARsolo/Solo.out/GeneFull/raw/matrix.mtx"
    output:
        GENEMAT = "{OUTDIR}/{SAMPLE}/rRNA/STARsolo/Solo.out/GeneFull/raw/matrix.mtx.gz"
    params:
        GENEDIR = directory("{OUTDIR}/{SAMPLE}/rRNA/STARsolo/Solo.out/GeneFull")
    threads:
        config['CORES']        
    run:
        recipe = RECIPE_DICT[wildcards.SAMPLE]
        if "noTrim" in recipe:
            #["seeker_v3.1_noTrimMatchLinker","seeker_v3.1_noTrim_total"]:
            shell(
                f"""
                cat {params.GENEDIR}/raw/barcodes.tsv | sed 's/_//' > {params.GENEDIR}/raw/barcodes_noUnderscore.tsv
                cat {params.GENEDIR}/filtered/barcodes.tsv | sed 's/_//' > {params.GENEDIR}/filtered/barcodes_noUnderscore.tsv
                """
            )

        # compress
        shell(
            f"""
            pigz -p{threads} {params.GENEDIR}/*/*.tsv {params.GENEDIR}/*/*.mtx 
            """
        )

# Index the .bam file produced by STAR
rule STAR_rRNA_indexSortedBAM:
    input:
        SORTEDBAM = '{OUTDIR}/{SAMPLE}/rRNA/STARsolo/Aligned.sortedByCoord.out.bam'
    output:
        BAI = '{OUTDIR}/{SAMPLE}/rRNA/STARsolo/Aligned.sortedByCoord.out.bam.bai'
    threads:
        config['CORES']
    run:
        shell(
            f"""
            {EXEC['SAMTOOLS']} index -@ {threads} {input.SORTEDBAM}
            """
        )


# Switch names because of STAR weirdness
rule STAR_rRNA_rename_compress_unmapped:
    input:
        UNMAPPED1 = '{OUTDIR}/{SAMPLE}/rRNA/STARsolo/Unmapped.out.mate1',
        UNMAPPED2 = '{OUTDIR}/{SAMPLE}/rRNA/STARsolo/Unmapped.out.mate2'
    output:
        FILTERED1_FQ = '{OUTDIR}/{SAMPLE}/rRNA/STARsolo/final_filtered_R1.fq.gz',
        FILTERED2_FQ = '{OUTDIR}/{SAMPLE}/rRNA/STARsolo/final_filtered_R2.fq.gz'
    params:
    threads:
        config['CORES']
    run:
        shell(
            f"""
            mv {input.UNMAPPED1} {output.FILTERED2_FQ.replace('.gz','')}
            mv {input.UNMAPPED2} {output.FILTERED1_FQ.replace('.gz','')}

            {EXEC['PIGZ']} \
                -p{threads} \
                {output.FILTERED2_FQ.replace('.gz','')} \
                {output.FILTERED1_FQ.replace('.gz','')}
            """
        )


#  Run fastqc on unmapped reads;
rule STAR_rRNA_filtered_fastqc:
    input:
        FQ  = '{OUTDIR}/{SAMPLE}/rRNA/STARsolo/final_filtered_{READ}.fq.gz',
    output:
        FQC_DIR = directory('{OUTDIR}/{SAMPLE}/fastqc/rRNA_STAR_{READ}')
    params:
        adapters = config['FASTQC_ADAPTERS']
    threads:
        config['CORES']
    run:
        shell(
            f"""
            mkdir -p {output.FQC_DIR}

            {EXEC['FASTQC']} \
                -o {output.FQC_DIR} \
                -t {threads} \
                -a {params.adapters} \
                {input.FQ}
            """
        )