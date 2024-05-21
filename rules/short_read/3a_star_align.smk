#############################################
## STAR Alignment
#############################################
# Make output directory, align fastqs, and generate raw/filtered feature/cell-barcode matrices
#   Info for STARsolo command line paramaters: https://github.com/alexdobin/STAR/blob/master/docs/STARsolo.md
# Run STARsolo
# TODO?: --twopassMode
# WASP?
rule STARsolo_align:
    input:
        FQS=lambda w: get_fqs(w, return_type="list", mode="ILMN"),
        WHITELIST=lambda w: get_whitelist(w, return_type="list"),
    output:
        BAM="{OUTDIR}/{SAMPLE}/STARsolo/short_read/{RECIPE}/Aligned.sortedByCoord.out.bam",  #TODO: add temp()
        UNMAPPED1="{OUTDIR}/{SAMPLE}/STARsolo/short_read/{RECIPE}/Unmapped.out.mate1",
        UNMAPPED2="{OUTDIR}/{SAMPLE}/STARsolo/short_read/{RECIPE}/Unmapped.out.mate2",
        VEL_MAT="{OUTDIR}/{SAMPLE}/STARsolo/short_read/{RECIPE}/Solo.out/Velocyto/raw/spliced.mtx",
        VEL_BC="{OUTDIR}/{SAMPLE}/STARsolo/short_read/{RECIPE}/Solo.out/Velocyto/raw/barcodes.tsv",
        VEL_FEAT="{OUTDIR}/{SAMPLE}/STARsolo/short_read/{RECIPE}/Solo.out/Velocyto/raw/features.tsv",
        GENE_MAT="{OUTDIR}/{SAMPLE}/STARsolo/short_read/{RECIPE}/Solo.out/Gene/raw/matrix.mtx",
        GENE_MAT_EM="{OUTDIR}/{SAMPLE}/STARsolo/short_read/{RECIPE}/Solo.out/Gene/raw/UniqueAndMult-EM.mtx",
        GENE_BC="{OUTDIR}/{SAMPLE}/STARsolo/short_read/{RECIPE}/Solo.out/Gene/raw/barcodes.tsv",
        GENE_FEAT="{OUTDIR}/{SAMPLE}/STARsolo/short_read/{RECIPE}/Solo.out/Gene/raw/features.tsv",
        GENEFULL_MAT="{OUTDIR}/{SAMPLE}/STARsolo/short_read/{RECIPE}/Solo.out/GeneFull/raw/matrix.mtx",
        GENEFULL_MAT_EM="{OUTDIR}/{SAMPLE}/STARsolo/short_read/{RECIPE}/Solo.out/GeneFull/raw/UniqueAndMult-EM.mtx",
        GENEFULL_BC="{OUTDIR}/{SAMPLE}/STARsolo/short_read/{RECIPE}/Solo.out/GeneFull/raw/barcodes.tsv",
        GENEFULL_FEAT="{OUTDIR}/{SAMPLE}/STARsolo/short_read/{RECIPE}/Solo.out/GeneFull/raw/features.tsv",
    params:
        STAR_REF=lambda w: get_STAR_ref(w),
        STAR_PARAMS=lambda w: get_STAR_extra_params(w),
    threads: config["CORES"]
    resources:
        MEMLIMIT=megabytes2bytes(config["MEMLIMIT_MB"]),
    priority: 42
    run:
        shell(
            f"""
            mkdir -p $(dirname {output.BAM})

            {EXEC['STAR']} \
                --runThreadN {threads} \
                --outFileNamePrefix $(dirname {output.BAM})/ \
                --outSAMtype BAM SortedByCoordinate \
                --limitBAMsortRAM={resources.MEMLIMIT} \
                --outSAMattributes NH HI nM AS CR UR CB UB GX GN sS sQ sM \
                --readFilesCommand zcat \
                --genomeDir {params.STAR_REF} \
                --readFilesIn {input.FQS[1]} {input.FQS[0]} \
                --outReadsUnmapped Fastx \
                --outSAMunmapped Within KeepPairs \
                --soloType {params.STAR_PARAMS['STAR.soloType']} {params.STAR_PARAMS['STAR.soloUMI']} {params.STAR_PARAMS['STAR.soloCB']} {params.STAR_PARAMS['STAR.soloAdapter']} {params.STAR_PARAMS['STAR.extra']} \
                --soloCBmatchWLtype {params.STAR_PARAMS['STAR.soloCBmatchWLtype']} \
                --soloCBwhitelist {" ".join(input.WHITELIST)} \
                --soloCellFilter TopCells $(wc -l {input.WHITELIST}) \
                --soloUMIfiltering MultiGeneUMI CR \
                --soloUMIdedup 1MM_CR \
                --soloBarcodeReadLength 0 \
                --soloFeatures Gene GeneFull Velocyto \
                --soloMultiMappers EM
            """
        )
                # --outSAMtype BAM SortedByCoordinate \
                # --clipAdapterType CellRanger4 \
        # --soloType {soloType} {soloUMI} {soloCB} {soloAdapter} {extraSTAR} \
        # --soloCBmatchWLtype {soloCBmatchWLtype} \
        # --soloCellFilter TopCells {nBB} \


# compress outputs from STAR (count matrices, cell barcodes, and gene lists)
rule compress_STAR_outs:
    input:
        VEL_MAT="{OUTDIR}/{SAMPLE}/STARsolo/short_read/{RECIPE}/Solo.out/Velocyto/raw/spliced.mtx",
        VEL_BC="{OUTDIR}/{SAMPLE}/STARsolo/short_read/{RECIPE}/Solo.out/Velocyto/raw/barcodes.tsv",
        VEL_FEAT="{OUTDIR}/{SAMPLE}/STARsolo/short_read/{RECIPE}/Solo.out/Velocyto/raw/features.tsv",
        GENE_MAT="{OUTDIR}/{SAMPLE}/STARsolo/short_read/{RECIPE}/Solo.out/Gene/raw/matrix.mtx",
        GENE_MAT_EM="{OUTDIR}/{SAMPLE}/STARsolo/short_read/{RECIPE}/Solo.out/Gene/raw/UniqueAndMult-EM.mtx",
        GENE_BC="{OUTDIR}/{SAMPLE}/STARsolo/short_read/{RECIPE}/Solo.out/Gene/raw/barcodes.tsv",
        GENE_FEAT="{OUTDIR}/{SAMPLE}/STARsolo/short_read/{RECIPE}/Solo.out/Gene/raw/features.tsv",
        GENEFULL_MAT="{OUTDIR}/{SAMPLE}/STARsolo/short_read/{RECIPE}/Solo.out/GeneFull/raw/matrix.mtx",
        GENEFULL_MAT_EM="{OUTDIR}/{SAMPLE}/STARsolo/short_read/{RECIPE}/Solo.out/GeneFull/raw/UniqueAndMult-EM.mtx",
        GENEFULL_BC="{OUTDIR}/{SAMPLE}/STARsolo/short_read/{RECIPE}/Solo.out/GeneFull/raw/barcodes.tsv",
        GENEFULL_FEAT="{OUTDIR}/{SAMPLE}/STARsolo/short_read/{RECIPE}/Solo.out/GeneFull/raw/features.tsv",
    output:
        VEL_MAT="{OUTDIR}/{SAMPLE}/STARsolo/short_read/{RECIPE}/Solo.out/Velocyto/raw/spliced.mtx.gz",
        VEL_BC="{OUTDIR}/{SAMPLE}/STARsolo/short_read/{RECIPE}/Solo.out/Velocyto/raw/barcodes.tsv.gz",
        VEL_FEAT="{OUTDIR}/{SAMPLE}/STARsolo/short_read/{RECIPE}/Solo.out/Velocyto/raw/features.tsv.gz",
        GENE_MAT="{OUTDIR}/{SAMPLE}/STARsolo/short_read/{RECIPE}/Solo.out/Gene/raw/matrix.mtx.gz",
        GENE_MAT_EM="{OUTDIR}/{SAMPLE}/STARsolo/short_read/{RECIPE}/Solo.out/Gene/raw/UniqueAndMult-EM.mtx.gz",
        GENE_BC="{OUTDIR}/{SAMPLE}/STARsolo/short_read/{RECIPE}/Solo.out/Gene/raw/barcodes.tsv.gz",
        GENE_FEAT="{OUTDIR}/{SAMPLE}/STARsolo/short_read/{RECIPE}/Solo.out/Gene/raw/features.tsv.gz",
        GENEFULL_MAT="{OUTDIR}/{SAMPLE}/STARsolo/short_read/{RECIPE}/Solo.out/GeneFull/raw/matrix.mtx.gz",
        GENEFULL_MAT_EM="{OUTDIR}/{SAMPLE}/STARsolo/short_read/{RECIPE}/Solo.out/GeneFull/raw/UniqueAndMult-EM.mtx.gz",
        GENEFULL_BC="{OUTDIR}/{SAMPLE}/STARsolo/short_read/{RECIPE}/Solo.out/GeneFull/raw/barcodes.tsv.gz",
        GENEFULL_FEAT="{OUTDIR}/{SAMPLE}/STARsolo/short_read/{RECIPE}/Solo.out/GeneFull/raw/features.tsv.gz",
    params:
        VELDIR=directory(
            "{OUTDIR}/{SAMPLE}/STARsolo/short_read/{RECIPE}/Solo.out/Velocyto"
        ),
        GENEDIR=directory(
            "{OUTDIR}/{SAMPLE}/STARsolo/short_read/{RECIPE}/Solo.out/Gene"
        ),
        GENEFULLDIR=directory(
            "{OUTDIR}/{SAMPLE}/STARsolo/short_read/{RECIPE}/Solo.out/GeneFull"
        ),
    threads: config["CORES"]
    run:
        if "noTrim" in wildcards.RECIPE and "seeker" in wildcards.RECIPE:
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
            {EXEC['PIGZ']} -p{threads} -f \
                {OUTDIR}/{wildcards.SAMPLE}/STARsolo/short_read/{wildcards.RECIPE}/*/*/*/*.tsv \
                {OUTDIR}/{wildcards.SAMPLE}/STARsolo/short_read/{wildcards.RECIPE}/*/*/*/*.mtx
            """
        )


# TODO
# rule multiqc_STAR:
#     input:
#         expand(
#             "{OUTDIR}/{SAMPLE}/STARsolo/short_read/{RECIPE}/Aligned.sortedByCoord.out.bam",
#             SAMPLE=SAMPLES, RECIPE=RECIPES),
#         expand("{OUTDIR}/{SAMPLE}/STARsolo/short_read/{RECIPE}/Unmapped.out.mate1",
#                SAMPLE=SAMPLES, RECIPE=RECIPES),
#         expand("{OUTDIR}/{SAMPLE}/STARsolo/short_read/{RECIPE}/Unmapped.out.mate2",
#                SAMPLE=SAMPLES, RECIPE=RECIPES),
#         expand("{OUTDIR}/{SAMPLE}/STARsolo/short_read/{RECIPE}/Solo.out/Velocyto/raw/spliced.mtx",
#                SAMPLE=SAMPLES, RECIPE=RECIPES),
#         expand("{OUTDIR}/{SAMPLE}/STARsolo/short_read/{RECIPE}/Solo.out/Velocyto/raw/barcodes.tsv",
#                SAMPLE=SAMPLES, RECIPE=RECIPES),
#         expand("{OUTDIR}/{SAMPLE}/STARsolo/short_read/{RECIPE}/Solo.out/Velocyto/raw/features.tsv",
#                SAMPLE=SAMPLES, RECIPE=RECIPES),
#         expand("{OUTDIR}/{SAMPLE}/STARsolo/short_read/{RECIPE}/Solo.out/Gene/raw/matrix.mtx",
#                SAMPLE=SAMPLES, RECIPE=RECIPES),
#         expand("{OUTDIR}/{SAMPLE}/STARsolo/short_read/{RECIPE}/Solo.out/Gene/raw/UniqueAndMult-EM.mtx",
#                SAMPLE=SAMPLES, RECIPE=RECIPES),
#         expand("{OUTDIR}/{SAMPLE}/STARsolo/short_read/{RECIPE}/Solo.out/Gene/raw/barcodes.tsv",
#                SAMPLE=SAMPLES, RECIPE=RECIPES),
#         expand("{OUTDIR}/{SAMPLE}/STARsolo/short_read/{RECIPE}/Solo.out/Gene/raw/features.tsv",
#                SAMPLE=SAMPLES, RECIPE=RECIPES),
#         expand("{OUTDIR}/{SAMPLE}/STARsolo/short_read/{RECIPE}/Solo.out/GeneFull/raw/matrix.mtx",
#                SAMPLE=SAMPLES, RECIPE=RECIPES),
#         expand("{OUTDIR}/{SAMPLE}/STARsolo/short_read/{RECIPE}/Solo.out/GeneFull/raw/UniqueAndMult-EM.mtx",
#                SAMPLE=SAMPLES, RECIPE=RECIPES),
#         expand("{OUTDIR}/{SAMPLE}/STARsolo/short_read/{RECIPE}/Solo.out/GeneFull/raw/barcodes.tsv",
#                SAMPLE=SAMPLES, RECIPE=RECIPES),
#         expand("{OUTDIR}/{SAMPLE}/STARsolo/short_read/{RECIPE}/Solo.out/GeneFull/raw/features.tsv",
#                SAMPLE=SAMPLES, RECIPE=RECIPES)
#     output:
#         REPORT = "{OUTDIR}/{SAMPLE}/STARsolo/short_read/multiqc_report.html"
#     params:
#         config['MEMLIMIT']
#     threads:
#         config['CORES']
#     resources:
#         mem_mb = config['MEMLIMIT_MB']
#     shell:
#         """
#         multiqc . -o {output.REPORT}
#         """
