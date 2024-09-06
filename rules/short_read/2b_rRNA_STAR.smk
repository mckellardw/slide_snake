# Make output directory, align fastqs, and generate raw/filtered feature/cell-barcode matrices
#   Info for STARsolo command line paramaters: https://github.com/alexdobin/STAR/blob/master/docs/STARsolo.md


# TODO- refactor to incorporate internal trimming options into rRNA filtering
## Need to run per (used) whitelist, per sample
rule ilmn_2b_STAR_rRNA_align:
    input:
        R1_FQ="{OUTDIR}/{SAMPLE}/tmp/twiceCut_R1.fq.gz",
        R2_FQ="{OUTDIR}/{SAMPLE}/tmp/twiceCut_R2.fq.gz",
        # WHITELIST=lambda w: get_whitelist(w, mode="all_used"),
        BC_WHITELIST="{OUTDIR}/{SAMPLE}/bc/whitelist.txt",
        BC_1="{OUTDIR}/{SAMPLE}/bc/whitelist_1.txt",
        BC_2="{OUTDIR}/{SAMPLE}/bc/whitelist_2.txt",
        BC_ADAPTER="{OUTDIR}/{SAMPLE}/bc/whitelist_adapter.txt",
        BC_US="{OUTDIR}/{SAMPLE}/bc/whitelist_underscore.txt",
    output:
        BAM="{OUTDIR}/{SAMPLE}/rRNA/STARsolo/Aligned.sortedByCoord.out.bam",  #TODO: add temp()
        UNMAPPED1="{OUTDIR}/{SAMPLE}/rRNA/STARsolo/Unmapped.out.mate1",
        UNMAPPED2="{OUTDIR}/{SAMPLE}/rRNA/STARsolo/Unmapped.out.mate2",
        GENEDIRECTORY=directory("{OUTDIR}/{SAMPLE}/rRNA/STARsolo/Solo.out/GeneFull"),
        GENEMAT="{OUTDIR}/{SAMPLE}/rRNA/STARsolo/Solo.out/GeneFull/raw/matrix.mtx",
    params:
        WHITELIST=lambda w: get_default_whitelist(w),
        STAR_REF=lambda w: get_STAR_ref(w, mode="rRNA"),
        STAR_PARAMS=lambda w: get_STAR_extra_params(w),
    resources:
        mem=megabytes2bytes(config["MEMLIMIT_MB"]),
        time="2:00:00",
    threads: config["CORES"]
    run:
        shell(
            f"""
            mkdir -p $(dirname {output.BAM})

            STAR \
                --runThreadN {threads} \
                --outFileNamePrefix $(dirname {output.BAM})/ \
                --outSAMtype BAM SortedByCoordinate \
                --outSAMattributes NH HI nM AS CR UR CB UB GX GN sS sQ sM \
                --readFilesCommand zcat \
                --genomeDir {params.STAR_REF} \
                --limitBAMsortRAM={resources.mem} \
                --readFilesIn {input.R2_FQ} {input.R1_FQ} \
                --outReadsUnmapped Fastx \
                --soloType {params.STAR_PARAMS['STAR.soloType']} {params.STAR_PARAMS['STAR.soloUMI']} {params.STAR_PARAMS['STAR.soloCB']} {params.STAR_PARAMS['STAR.soloAdapter']} {params.STAR_PARAMS['STAR.extra']} \
                --soloCBwhitelist {params.WHITELIST} \
                --soloCellFilter TopCells $(wc -l {params.WHITELIST}) \
                --soloCBmatchWLtype {params.STAR_PARAMS['STAR.soloCBmatchWLtype']} \
                --soloUMIfiltering MultiGeneUMI CR \
                --soloUMIdedup 1MM_CR \
                --soloBarcodeReadLength 0 \
                --soloFeatures GeneFull \
                --soloMultiMappers EM
            """
        )
        # --outSAMunmapped Within KeepPairs \
        # --clipAdapterType CellRanger4 \



# compress outputs from STAR (count matrices, cell barcodes, and gene lists)
rule ilmn_2b_STAR_rRNA_compress_outs:
    input:
        GENEMAT="{OUTDIR}/{SAMPLE}/rRNA/STARsolo/Solo.out/GeneFull/raw/matrix.mtx",
    output:
        GENEMAT="{OUTDIR}/{SAMPLE}/rRNA/STARsolo/Solo.out/GeneFull/raw/matrix.mtx.gz",
    params:
        GENEDIR=directory("{OUTDIR}/{SAMPLE}/rRNA/STARsolo/Solo.out/GeneFull"),
    threads: config["CORES"]
    run:
        recipe = RECIPE_DICT[wildcards.SAMPLE]
        if "noTrim" in recipe:
            # ["seeker_noTrimMatchLinker","seeker_noTrim_total"]:
            shell(
                f"""
                cat {params.GENEDIR}/raw/barcodes.tsv \
                | sed 's/_//' \
                > {params.GENEDIR}/raw/barcodes_noUnderscore.tsv
                
                cat {params.GENEDIR}/filtered/barcodes.tsv \
                | sed 's/_//' \
                > {params.GENEDIR}/filtered/barcodes_noUnderscore.tsv
                """
            )

            # compress
        shell(
            f"""
            pigz \
                -p{threads} \
                {params.GENEDIR}/*/*.tsv \
                {params.GENEDIR}/*/*.mtx 
            """
        )


# Switch names because of STAR weirdness
rule ilmn_2b_STAR_rRNA_rename_compress_unmapped:
    input:
        UNMAPPED1="{OUTDIR}/{SAMPLE}/rRNA/STARsolo/Unmapped.out.mate1",
        UNMAPPED2="{OUTDIR}/{SAMPLE}/rRNA/STARsolo/Unmapped.out.mate2",
    output:
        FILTERED1_FQ="{OUTDIR}/{SAMPLE}/rRNA/STARsolo/noRibo_R1.fq.gz",
        FILTERED2_FQ="{OUTDIR}/{SAMPLE}/rRNA/STARsolo/noRibo_R2.fq.gz",
    threads: config["CORES"]
    run:
        shell(
            f"""
            mv {input.UNMAPPED1} {output.FILTERED2_FQ.replace('.gz' , '')}
            mv {input.UNMAPPED2} {output.FILTERED1_FQ.replace('.gz' , '')}

            pigz \
                -p{threads} \
                {output.FILTERED2_FQ.replace('.gz' , '')} \
                {output.FILTERED1_FQ.replace('.gz' , '')}
            """
        )


#  Run fastqc on unmapped reads;
rule ilmn_2b_STAR_rRNA_filtered_fastqc:
    input:
        FQ="{OUTDIR}/{SAMPLE}/rRNA/STARsolo/noRibo_{READ}.fq.gz",
    output:
        FQC_DIR=directory("{OUTDIR}/{SAMPLE}/fastqc/rRNA_STAR_{READ}"),
    params:
        adapters=config["FASTQC_ADAPTERS"],
    threads: config["CORES"]
    conda:
        f"{workflow.basedir}/envs/fastqc.yml"
    shell:
        """
        mkdir -p {output.FQC_DIR}

        fastqc \
            -o {output.FQC_DIR} \
            -t {threads} \
            -a {params.adapters} \
            {input.FQ}
        """
