# Filter out rRNA reads w/ bwa alignment
# VASAseq implementation - https://github.com/anna-alemany/VASAseq/blob/main/mapping/ribo-bwamem.sh
rule bwa_rRNA_align:
    input:
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
        # MEMLIMIT = config['MEMLIMIT'],
        BWA_REF = lambda w: rRNA_BWA_DICT[w.SAMPLE],
        MIN_ALIGNSCORE = 40,
    log:
        log = "{OUTDIR}/{SAMPLE}/rRNA/bwa/bwa_mem.log"
    threads:
        config['CORES']
    run:
        # Align to rRNA ref w/ `bwa mem` for cleaner/faster rRNA filtering 
        ##TODO incorporate VASAseq style "long"/short read handling with multiple align steps
        ##TODO- modify this so that it doesn't treat reads as paired-end
        shell(            
            f"""
            mkdir -p $(dirname {output.BAM1})

            {EXEC['BWA']} mem \
                -t {threads} \
                {params.BWA_REF} \
                {input.R2_FQ} \
            1> {output.BAM1} \
            2> {log.log} \
            
            {EXEC['SAMTOOLS']} sort \
                -@ {threads} \
                -O BAM \
                {output.BAM1} \
            > {output.BAM2} 

            {EXEC['SAMTOOLS']} view \
                -h {output.BAM2} \
            | awk \
                -v quality={params.MIN_ALIGNSCORE} \
                -f scripts/awk/bam_lowPassMAPQFilter.awk \
            | {EXEC['SAMTOOLS']} fastq \
                {output.BAM2} \
            > {output.R2_FQ_BWA_FILTERED} 
            """
        )


rule bwa_rRNA_index_alignment:
    input:
        BAM = "{OUTDIR}/{SAMPLE}/rRNA/bwa/aligned_sorted.bam", 
    output:
        BAI = "{OUTDIR}/{SAMPLE}/rRNA/bwa/aligned_sorted.bam.bai",
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
            > {output.R1_FQ_BWA_FILTERED}
            """
        )


rule bwa_rRNA_compress_unmapped:
    input:
        R1_FQ_BWA_FILTERED  = '{OUTDIR}/{SAMPLE}/rRNA/bwa/final_filtered_R1.fq',
        R2_FQ_BWA_FILTERED  = '{OUTDIR}/{SAMPLE}/rRNA/bwa/final_filtered_R2.fq'
    output:        
        R1_FQ_BWA_FILTERED  = '{OUTDIR}/{SAMPLE}/rRNA/bwa/final_filtered_R1.fq.gz',
        R2_FQ_BWA_FILTERED  = '{OUTDIR}/{SAMPLE}/rRNA/bwa/final_filtered_R2.fq.gz'
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