# Make output directory, align fastqs, and generate raw/filtered feature/cell-barcode matrices
#   Info for STARsolo command line paramaters: https://github.com/alexdobin/STAR/blob/master/docs/STARsolo.md
#TODO- dedup rRNA .bam files
rule bwa_align_rRNA:
    input:
        uBAM = temp('{OUTDIR}/{SAMPLE}/tmp/unaligned_barcoded.bam')
        BB_WHITELIST = "{OUTDIR}/{SAMPLE}/bb/whitelist.txt",
        BB_1 = "{OUTDIR}/{SAMPLE}/bb/whitelist_1.txt",
        BB_2 = "{OUTDIR}/{SAMPLE}/bb/whitelist_2.txt",
        BB_ADAPTER = "{OUTDIR}/{SAMPLE}/bb/whitelist_adapter.txt"
    output:
        UNMAPPED1 = '{OUTDIR}/{SAMPLE}/rRNA/bwa/unmapped_R1.fq',
        UNMAPPED2 = '{OUTDIR}/{SAMPLE}/rRNA/bwa/unmapped_R2.fq'
    params:
        MEMLIMIT = config['MEMLIMIT']
    log:
        log = "{OUTDIR}/{SAMPLE}/rRNA/bwa_rRNA.log"
    threads:
        config['CORES']
    run: 
        # tmp_recipe = RECIPE_DICT[wildcards.SAMPLE]
        # STAR_REF = rRNA_DICT[wildcards.SAMPLE] # use rRNA ref
        # nBB = sum(1 for line in open(input.BB_WHITELIST)) # get number of bead barcodes for filtered count matrix, `--soloCellFilter`

        # soloType = RECIPE_SHEET["STAR.soloType"][tmp_recipe]
        # soloUMI = RECIPE_SHEET["STAR.soloUMI"][tmp_recipe]
        # soloCB = RECIPE_SHEET["STAR.soloCB"][tmp_recipe]
        # soloCBmatchWLtype = RECIPE_SHEET["STAR.soloCBmatchWLtype"][tmp_recipe]
        # soloAdapter = RECIPE_SHEET["STAR.soloAdapter"][tmp_recipe]
        # extraSTAR = RECIPE_SHEET["STAR.extra"][tmp_recipe]

        if "bwa" in tmp_recipe:
            # Align to rRNA ref w/ `bwa mem` for cleaner/faster rRNA filtering 
            ## Add empty count mat file for 
            ##TODO: fix .fasta passing here (*.fa is dangerous...)
            ##TODO incorporate VASAseq style "long"/short read handling with multiple align steps
            shell(
                f"""
                {EXEC['BWA']} mem \
                    -t {threads} \
                    {STAR_REF}/*.fa \
                    {input.R1_FQ} {input.R2_FQ} \
                | {EXEC['SAMTOOLS']} sort \
                    -@ {threads} \
                | {EXEC['SAMTOOLS']} fastq \
                    -f 4 \
                    -1 {output.UNMAPPED2} \
                    -2 {output.UNMAPPED1} \
                    -0 /dev/null \
                | tee {log.log}
                """
            )


                # {EXEC['BWA']} mem \
                #     -t {threads} \
                #     {STAR_REF}/*.fa \
                #     {input.R1_FQ} {input.R2_FQ} \
                # > {OUTDIR}/{wildcards.SAMPLE}/STARsolo_rRNA/Aligned.sortedByCoord.out.sam \
                # | tee {log.log}

                # {EXEC['SAMTOOLS']} sort \
                #     -@ {threads} \
                #     {OUTDIR}/{wildcards.SAMPLE}/STARsolo_rRNA/Aligned.sortedByCoord.out.sam \
                # > {OUTDIR}/{wildcards.SAMPLE}/STARsolo_rRNA/Aligned.sortedByCoord.out.bam

                # {EXEC['SAMTOOLS']} fastq \
                #     -f 4 -1 {output.UNMAPPED2} \
                #     -2 {output.UNMAPPED1} \
                #     -0 /dev/null {OUTDIR}/{wildcards.SAMPLE}/STARsolo_rRNA/Aligned.sortedByCoord.out.bam


rule compress_unmapped_bwa:
    input:
        UNMAPPED1 = '{OUTDIR}/{SAMPLE}/rRNA/bwa/unmapped_R1.fq',
        UNMAPPED2 = '{OUTDIR}/{SAMPLE}/rRNA/bwa/unmapped_R2.fq'
    output:
        UNMAPPED1 = '{OUTDIR}/{SAMPLE}/rRNA/bwa/unmapped_R1.fq.gz',
        UNMAPPED2 = '{OUTDIR}/{SAMPLE}/rRNA/bwa/unmapped_R2.fq.gz'
    params:
    threads:
        config['CORES']
    run:
        shell(
            f"""
            {EXEC['PIGZ']} -p{threads} {input}
            """
        )


#  Run fastqc on unmapped reads;
rule rRNA_filtered_fastqc_bwa:
    input:
        FILTERED1_FQ = '{OUTDIR}/{SAMPLE}/rRNA/bwa/unmapped_R1.fq.gz',
        FILTERED2_FQ = '{OUTDIR}/{SAMPLE}/rRNA/bwa/unmapped_R2.fq.gz'
    output:
        FQC_DIR = directory('{OUTDIR}/{SAMPLE}/fastqc/rRNA_filtered_bwa')
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
                {input.FILTERED1_FQ} {input.FILTERED2_FQ}
            """
        )