# Make output directory, align fastqs, and generate raw/filtered feature/cell-barcode matrices
#   Info for STARsolo command line paramaters: https://github.com/alexdobin/STAR/blob/master/docs/STARsolo.md

rule STARsolo_align_rRNA:
    input:
        R1_FQ_HardTrim = '{OUTDIR}/{sample}/tmp/{sample}_R1_final.fq.gz',
        R1_FQ = '{OUTDIR}/{sample}/tmp/{sample}_R1_final.fq.gz',
        R2_FQ = '{OUTDIR}/{sample}/tmp/{sample}_R2_final.fq.gz',
        BB_WHITELIST = "{OUTDIR}/{sample}/bb/whitelist.txt",
        BB_1 = "{OUTDIR}/{sample}/bb/whitelist_1.txt",
        BB_2 = "{OUTDIR}/{sample}/bb/whitelist_2.txt"
    output:
        SORTEDBAM = '{OUTDIR}/{sample}/STARsolo_rRNA/Aligned.sortedByCoord.out.bam', #TODO: add temp()
        UNMAPPED1 = '{OUTDIR}/{sample}/STARsolo_rRNA/Unmapped.out.mate1',
        UNMAPPED2 = '{OUTDIR}/{sample}/STARsolo_rRNA/Unmapped.out.mate2',
        GENE = directory('{OUTDIR}/{sample}/STARsolo_rRNA/Solo.out/GeneFull'),
        GENEMAT = '{OUTDIR}/{sample}/STARsolo_rRNA/Solo.out/GeneFull/raw/matrix.mtx'
    params:
        STAR_EXEC = config['STAR_EXEC'],
        MEMLIMIT = config['MEMLIMIT']
    threads:
        config['CORES']
    run: 
        tmp_chemistry = CHEM_DICT[wildcards.sample]
        STAR_REF = rRNA_DICT[wildcards.sample] # use rRNA ref
        
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

        R2 = input.R2_FQ

        shell(
            f"""
            mkdir -p {OUTDIR}/{wildcards.sample}/STARsolo_rRNA

            {params.STAR_EXEC} \
            --runThreadN {threads} \
            --outFileNamePrefix {OUTDIR}/{wildcards.sample}/STARsolo_rRNA/ \
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
            --soloFeatures GeneFull \
            --soloMultiMappers EM
            """
        )

# compress outputs from STAR (count matrices, cell barcodes, and gene lists)
rule compress_STAR_rRNA_outs:
    input:
        GENEMAT = "{OUTDIR}/{sample}/STARsolo_rRNA/Solo.out/GeneFull/raw/matrix.mtx"
    output:
        GENEMAT = "{OUTDIR}/{sample}/STARsolo_rRNA/Solo.out/GeneFull/raw/matrix.mtx.gz"
    params:
        GENEDIR = directory("{OUTDIR}/{sample}/STARsolo_rRNA/Solo.out/GeneFull")
    threads:
        config['CORES']        
    run:
        shell(
            f"""
            pigz -p{threads} {params.GENEDIR}/*/*.tsv {params.GENEDIR}/*/*.mtx 
            """
        )

# Index the .bam file produced by STAR
rule indexSortedBAM_rRNA:
    input:
        SORTEDBAM = '{OUTDIR}/{sample}/STARsolo_rRNA/Aligned.sortedByCoord.out.bam'
    output:
        BAI = '{OUTDIR}/{sample}/STARsolo_rRNA/Aligned.sortedByCoord.out.bam.bai'
    threads:
        config['CORES']
    run:
        shell(
            f"""
            samtools index -@ {threads} {input.SORTEDBAM}
            """
        )


# Run fastqc on unmapped reads; switch names because of STAR weirdness
rule rRNA_filtered_fastqc:
    input:
        UNMAPPED1 = '{OUTDIR}/{sample}/STARsolo_rRNA/Unmapped.out.mate1',
        UNMAPPED2 = '{OUTDIR}/{sample}/STARsolo_rRNA/Unmapped.out.mate2'
    output:
        FILTERED1_FQ = '{OUTDIR}/{sample}/tmp/{sample}_R1_final_filtered.fq.gz',
        FILTERED2_FQ = '{OUTDIR}/{sample}/tmp/{sample}_R2_final_filtered.fq.gz',
        FQC_DIR = directory('{OUTDIR}/{sample}/rRNA_filtered_fastqc_out')
    params:
        FASTQC_EXEC = config['FASTQC_EXEC']
    threads:
        config['CORES']
    run:
        shell(
            f"""
            mv {input.UNMAPPED1} {OUTDIR}/{wildcards.sample}/tmp/{wildcards.sample}_R2_final_filtered.fq
            mv {input.UNMAPPED2} {OUTDIR}/{wildcards.sample}/tmp/{wildcards.sample}_R1_final_filtered.fq

            pigz -p{threads} {OUTDIR}/{wildcards.sample}/tmp/*.fq

            mkdir -p {output.FQC_DIR}

            {params.FASTQC_EXEC} \
             -o {output.FQC_DIR} \
             -t {threads} \
             {output.FILTERED1_FQ} {output.FILTERED2_FQ}
            """
        )