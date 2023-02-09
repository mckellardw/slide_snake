# Make output directory, align fastqs, and generate raw/filtered feature/cell-barcode matrices
#   Info for STARsolo command line paramaters: https://github.com/alexdobin/STAR/blob/master/docs/STARsolo.md

rule STARsolo_align_rRNA:
    input:
        FINAL_R1_FQ = '{OUTDIR}/{sample}/tmp/{sample}_R1_final.fq.gz',
        FINAL_R2_FQ = '{OUTDIR}/{sample}/tmp/{sample}_R2_final.fq.gz'
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

        UMIlen = CHEMISTRY_SHEET["STAR.UMIlen"][tmp_chemistry]
        SOLOtype = CHEMISTRY_SHEET["STAR.soloType"][tmp_chemistry]
        CB_WHITELIST = CHEMISTRY_SHEET["whitelist"][tmp_chemistry]
        # extraSTAR = CHEMISTRY_SHEET["STAR.extra"][tmp_chemistry]

        shell(
            f"""
            mkdir -p {OUTDIR}/{wildcards.sample}

            {params.STAR_EXEC} \
            --runThreadN {threads} \
            --outFileNamePrefix {OUTDIR}/{wildcards.sample}/STARsolo_rRNA/ \
            --outSAMtype BAM SortedByCoordinate \
            --outSAMattributes NH HI nM AS CR UR CB UB GX GN sS sQ sM \
            --readFilesCommand zcat \
            --soloUMIlen {UMIlen} \
            --genomeDir {STAR_REF} \
            --genomeLoad LoadAndRemove \
            --limitBAMsortRAM={params.MEMLIMIT} \
            --readFilesIn {input.FINAL_R2_FQ} {input.FINAL_R1_FQ} \
            --clipAdapterType CellRanger4 \
            --outReadsUnmapped Fastx \
            --soloType {SOLOtype} \
            --soloBarcodeReadLength 0 \
            --soloCBwhitelist {CB_WHITELIST} \
            --soloCellFilter EmptyDrops_CR \
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