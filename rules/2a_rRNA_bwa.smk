# Filter out rRNA reads w/ bwa alignment
# VASAseq implementation - https://github.com/anna-alemany/VASAseq/blob/main/mapping/ribo-bwamem.sh
rule bwa_rRNA_align:
    input:
        # uBAM = temp('{OUTDIR}/{SAMPLE}/tmp/unaligned_barcoded.bam'),
        # R2_FQ = '{OUTDIR}/{SAMPLE}/tmp/cut_R2.fq.gz',
        R2_FQ = "{OUTDIR}/{SAMPLE}/tmp/twiceCut_R2.fq.gz",
        BB_WHITELIST = "{OUTDIR}/{SAMPLE}/bb/whitelist.txt",
        BB_1 = "{OUTDIR}/{SAMPLE}/bb/whitelist_1.txt",
        BB_2 = "{OUTDIR}/{SAMPLE}/bb/whitelist_2.txt",
        BB_ADAPTER = "{OUTDIR}/{SAMPLE}/bb/whitelist_adapter.txt"
    output:
        BAM1 = temp("{OUTDIR}/{SAMPLE}/rRNA/bwa/aligned.bam"),
        BAM2 = "{OUTDIR}/{SAMPLE}/rRNA/bwa/aligned_sorted.bam", 
        # R1_FQ_BWA_FILTERED  = '{OUTDIR}/{SAMPLE}/rRNA/bwa/final_filtered_R1.fq',
        R2_FQ_BWA_FILTERED  = "{OUTDIR}/{SAMPLE}/rRNA/bwa/final_filtered_R2.fq"
    params:
        MEMLIMIT = config['MEMLIMIT']
    log:
        log = "{OUTDIR}/{SAMPLE}/rRNA/bwa/bwa_rRNA.log"
    threads:
        config['CORES']
    run: 
        # tmp_recipe = RECIPE_DICT[wildcards.SAMPLE]
        BWA_REF = rRNA_BWA_DICT[wildcards.SAMPLE] # use rRNA ref
        # nBB = sum(1 for line in open(input.BB_WHITELIST)) # get number of bead barcodes for filtered count matrix, `--soloCellFilter`

        # soloType = RECIPE_SHEET["STAR.soloType"][tmp_recipe]
        # soloUMI = RECIPE_SHEET["STAR.soloUMI"][tmp_recipe]
        # soloCB = RECIPE_SHEET["STAR.soloCB"][tmp_recipe]
        # soloCBmatchWLtype = RECIPE_SHEET["STAR.soloCBmatchWLtype"][tmp_recipe]
        # soloAdapter = RECIPE_SHEET["STAR.soloAdapter"][tmp_recipe]
        # extraSTAR = RECIPE_SHEET["STAR.extra"][tmp_recipe]

        # Align to rRNA ref w/ `bwa mem` for cleaner/faster rRNA filtering 
        ##TODO incorporate VASAseq style "long"/short read handling with multiple align steps
        ##TODO- modify this so that it doesn't treat reads as paired-end
        shell(
            f"""
            mkdir -p $(dirname {output.BAM1})

            {EXEC['BWA']} mem \
                -t {threads} \
                {BWA_REF} \
                {input.R2_FQ} \
            1> {output.BAM1} \
            2> {log.log} \
            
            {EXEC['SAMTOOLS']} sort \
                -@ {threads} \
                -O BAM \
                {output.BAM1} \
            > {output.BAM2} 
            
            {EXEC['SAMTOOLS']} fastq \
                -f 4 \
                {output.BAM2} \
            > {output.R2_FQ_BWA_FILTERED} 
            """
        )

                # -2 \
                # -0 /dev/null \
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

rule bwa_rRNA_index_alignment:
    input:
        BAM = '{OUTDIR}/{SAMPLE}/rRNA/bwa/aligned_sorted.bam', 
    output:
        BAI = '{OUTDIR}/{SAMPLE}/rRNA/bwa/aligned_sorted.bam.bai',
    run:
        shell(
            f"""
            {EXEC['SAMTOOLS']} index {input.BAM}
            """
        )


rule bwa_rRNA_filter_R1:
    input:
        R1_FQ = '{OUTDIR}/{SAMPLE}/tmp/cut_R1.fq.gz',
        R2_FQ = '{OUTDIR}/{SAMPLE}/tmp/cut_R2.fq.gz',
        R2_FQ_BWA_FILTERED  = '{OUTDIR}/{SAMPLE}/rRNA/bwa/final_filtered_R2.fq'
    output:
        R1_FQ = temp('{OUTDIR}/{SAMPLE}/tmp/cut_R1.fq'),
        R1_FQ_BWA_FILTERED  = '{OUTDIR}/{SAMPLE}/rRNA/bwa/final_filtered_R1.fq',
        # rRNA_LIST  = '{OUTDIR}/{SAMPLE}/rRNA/bwa/rRNA_readID.list', #temp()
        UNMAPPED_LIST  = '{OUTDIR}/{SAMPLE}/rRNA/bwa/rRNA_readID.list'  #temp()
    params:
    threads:
        config['CORES']
    run:
        shell(
            f"""
            cat {input.R2_FQ_BWA_FILTERED} \
            | awk -f scripts/awk/fq_readHeader.awk - \
            > {output.UNMAPPED_LIST}

            zcat {input.R1_FQ} > {output.R1_FQ} 
            
            {EXEC['SEQTK']} subseq {output.R1_FQ} {output.UNMAPPED_LIST} \
            > {output.R1_FQ_BWA_FILTERED.strip('.gz')}
            """
        )

        # {EXEC['PIGZ']} -p{threads} {output.R1_FQ_BWA_FILTERED.strip('.gz')}


rule bwa_rRNA_compress_unmapped:
    input:
        R1_FQ_BWA_FILTERED  = '{OUTDIR}/{SAMPLE}/rRNA/bwa/final_filtered_R1.fq',
        R2_FQ_BWA_FILTERED  = '{OUTDIR}/{SAMPLE}/rRNA/bwa/final_filtered_R2.fq'
    output:        
        R1_FQ_BWA_FILTERED  = '{OUTDIR}/{SAMPLE}/rRNA/bwa/final_filtered_R1.fq.gz',
        R2_FQ_BWA_FILTERED  = '{OUTDIR}/{SAMPLE}/rRNA/bwa/final_filtered_R2.fq.gz'
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
rule bwa_rRNA_filtered_fastqc:
    input:
        FQ_BWA_FILTERED  = '{OUTDIR}/{SAMPLE}/rRNA/bwa/final_filtered_{READ}.fq.gz',
    output:
        FQC_DIR = directory('{OUTDIR}/{SAMPLE}/fastqc/rRNA_bwa_{READ}')
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
                {input.FQ_BWA_FILTERED}
            """
        )