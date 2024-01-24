# Make output directory, align fastqs, and generate raw/filtered feature/cell-barcode matrices
#   Info for STARsolo command line paramaters: https://github.com/alexdobin/STAR/blob/master/docs/STARsolo.md
#TODO- dedup rRNA .bam files
rule STARsolo_align_rRNA:
    input:
        R1_FQ = '{OUTDIR}/{sample}/tmp/final_R1.fq.gz',
        R2_FQ = '{OUTDIR}/{sample}/tmp/final_R2.fq.gz',
        BB_WHITELIST = "{OUTDIR}/{sample}/bb/whitelist.txt",
        BB_1 = "{OUTDIR}/{sample}/bb/whitelist_1.txt",
        BB_2 = "{OUTDIR}/{sample}/bb/whitelist_2.txt",
        BB_ADAPTER = "{OUTDIR}/{sample}/bb/whitelist_adapter.txt"
    output:
        SORTEDBAM = '{OUTDIR}/{sample}/STARsolo_rRNA/Aligned.sortedByCoord.out.bam', #TODO: add temp()
        UNMAPPED1 = '{OUTDIR}/{sample}/STARsolo_rRNA/Unmapped.out.mate1',
        UNMAPPED2 = '{OUTDIR}/{sample}/STARsolo_rRNA/Unmapped.out.mate2',
        GENEDIRECTORY = directory('{OUTDIR}/{sample}/STARsolo_rRNA/Solo.out/GeneFull'),
        GENEMAT = '{OUTDIR}/{sample}/STARsolo_rRNA/Solo.out/GeneFull/raw/matrix.mtx'
    params:
        MEMLIMIT = config['MEMLIMIT']
    threads:
        config['CORES']
    run: 
        tmp_recipe = RECIPE_DICT[wildcards.sample]
        STAR_REF = rRNA_DICT[wildcards.sample] # use rRNA ref
        nBB = sum(1 for line in open(input.BB_WHITELIST)) # get number of bead barcodes for filtered count matrix, `--soloCellFilter`

        #TODO: add try catches
        soloType = RECIPE_SHEET["STAR.soloType"][tmp_recipe]
        soloUMI = RECIPE_SHEET["STAR.soloUMI"][tmp_recipe]
        soloCB = RECIPE_SHEET["STAR.soloCB"][tmp_recipe]
        soloCBmatchWLtype = RECIPE_SHEET["STAR.soloCBmatchWLtype"][tmp_recipe]
        soloAdapter = RECIPE_SHEET["STAR.soloAdapter"][tmp_recipe]
        extraSTAR = RECIPE_SHEET["STAR.extra"][tmp_recipe]

        # Params for different barcode handling strategies
        if "stomics" in tmp_recipe:
            whitelist = input.BB_WHITELIST
        elif "noTrim" in tmp_recipe:
            whitelist = f"{input.BB_1} {input.BB_2}"
        #     R1 = input.R1_FQ
        elif "internalTrim" in tmp_recipe:
            whitelist = input.BB_WHITELIST
        #     R1 = input.R1_FQ_InternalTrim
        elif "adapterInsert" in tmp_recipe:
            whitelist = input.BB_ADAPTER
        #     R1 = input.R1_FQ
        else:
            whitelist = input.BB_WHITELIST
        #     R1 = input.R1_FQ_HardTrim

        # R2 = input.R2_FQ

        if "bwa" in tmp_recipe:
            # Align to rRNA ref w/ `bwa mem` for cleaner/faster rRNA filtering 
            ## Flip read numbers to match STAR
            ## Add empty count mat file for 
            ##TODO: fix .fasta passing here (*.fa is dangerous...)
            ##TODO: fix folder/file naming schema...
            ##TODO: add log file
            shell(
                f"""
                {EXEC['BWA']} mem -t {threads} {STAR_REF}/*.fa {input.R1_FQ} {input.R2_FQ} \
                > {OUTDIR}/{wildcards.sample}/STARsolo_rRNA/Aligned.sortedByCoord.out.sam \
                | tee {OUTDIR}/{wildcards.sample}/STARsolo_rRNA/bwa.log

                {EXEC['SAMTOOLS']} sort -@ {threads} {OUTDIR}/{wildcards.sample}/STARsolo_rRNA/Aligned.sortedByCoord.out.sam \
                > {OUTDIR}/{wildcards.sample}/STARsolo_rRNA/Aligned.sortedByCoord.out.bam

                {EXEC['SAMTOOLS']} fastq -f 4 -1 {output.UNMAPPED2} -2 {output.UNMAPPED1} -0 /dev/null {OUTDIR}/{wildcards.sample}/STARsolo_rRNA/Aligned.sortedByCoord.out.bam

                mkdir -p {output.GENEDIRECTORY}
                echo "empty placeholder file..." > {output.GENEMAT} 
                """
            )
                # rm {OUTDIR}/{wildcards.sample}/STARsolo_rRNA/STARsolo_rRNA/Aligned.sortedByCoord.out.sam
                # | {EXEC['SAMTOOLS']} view -f 4 -bh - \
        else:
            shell(
                f"""
                mkdir -p {OUTDIR}/{wildcards.sample}/STARsolo_rRNA

                {EXEC['STAR']} \
                    --runThreadN {threads} \
                    --outFileNamePrefix {OUTDIR}/{wildcards.sample}/STARsolo_rRNA/ \
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
        tmp_recipe = RECIPE_DICT[wildcards.sample]
        if "noTrim" in tmp_recipe:
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
            {EXEC['SAMTOOLS']} index -@ {threads} {input.SORTEDBAM}
            """
        )


# Run fastqc on unmapped reads; switch names because of STAR weirdness
##TODO: split fastqc into separate rule
rule rRNA_filtered_fastqc:
    input:
        UNMAPPED1 = '{OUTDIR}/{sample}/STARsolo_rRNA/Unmapped.out.mate1',
        UNMAPPED2 = '{OUTDIR}/{sample}/STARsolo_rRNA/Unmapped.out.mate2'
    output:
        FILTERED1_FQ = '{OUTDIR}/{sample}/tmp/final_filtered_R1.fq.gz',
        FILTERED2_FQ = '{OUTDIR}/{sample}/tmp/final_filtered_R2.fq.gz',
        FQC_DIR = directory('{OUTDIR}/{sample}/fastqc/rRNA_filtered')
    params:
    threads:
        config['CORES']
    run:
        shell(
            f"""
            mv {input.UNMAPPED1} {OUTDIR}/{wildcards.sample}/tmp/{wildcards.sample}_R2_final_filtered.fq
            mv {input.UNMAPPED2} {OUTDIR}/{wildcards.sample}/tmp/{wildcards.sample}_R1_final_filtered.fq

            pigz -p{threads} -f {OUTDIR}/{wildcards.sample}/tmp/*.fq

            mkdir -p {output.FQC_DIR}

            {EXEC['FASTQC']} \
             -o {output.FQC_DIR} \
             -t {threads} \
             {output.FILTERED1_FQ} {output.FILTERED2_FQ}
            """
        )