#############################################
## STAR Alignment
#############################################
# Make output directory, align fastqs, and generate raw/filtered feature/cell-barcode matrices
#   Info for STARsolo command line paramaters: https://github.com/alexdobin/STAR/blob/master/docs/STARsolo.md
# Run STARsolo
# WASP?
# rule ilmn_3a_STARsolo_align:
#     input:
#         FQS=lambda w: get_fqs(w, return_type="list", mode="ILMN"),
#         WHITELIST=lambda w: get_whitelist(w, return_type="list"),
#     output:
#         BAM="{OUTDIR}/{SAMPLE}/short_read/STARsolo/{RECIPE}/Aligned.sortedByCoord.out.bam",  #TODO: add temp()
#         UNMAPPED1="{OUTDIR}/{SAMPLE}/short_read/STARsolo/{RECIPE}/Unmapped.out.mate1",
#         UNMAPPED2="{OUTDIR}/{SAMPLE}/short_read/STARsolo/{RECIPE}/Unmapped.out.mate2",
#         VEL_MAT="{OUTDIR}/{SAMPLE}/short_read/STARsolo/{RECIPE}/Solo.out/Velocyto/raw/spliced.mtx",
#         VEL_BC="{OUTDIR}/{SAMPLE}/short_read/STARsolo/{RECIPE}/Solo.out/Velocyto/raw/barcodes.tsv",
#         VEL_FEAT="{OUTDIR}/{SAMPLE}/short_read/STARsolo/{RECIPE}/Solo.out/Velocyto/raw/features.tsv",
#         GENE_MAT="{OUTDIR}/{SAMPLE}/short_read/STARsolo/{RECIPE}/Solo.out/Gene/raw/matrix.mtx",
#         GENE_MAT_EM="{OUTDIR}/{SAMPLE}/short_read/STARsolo/{RECIPE}/Solo.out/Gene/raw/UniqueAndMult-EM.mtx",
#         GENE_BC="{OUTDIR}/{SAMPLE}/short_read/STARsolo/{RECIPE}/Solo.out/Gene/raw/barcodes.tsv",
#         GENE_FEAT="{OUTDIR}/{SAMPLE}/short_read/STARsolo/{RECIPE}/Solo.out/Gene/raw/features.tsv",
#         GENEFULL_MAT="{OUTDIR}/{SAMPLE}/short_read/STARsolo/{RECIPE}/Solo.out/GeneFull/raw/matrix.mtx",
#         GENEFULL_MAT_EM="{OUTDIR}/{SAMPLE}/short_read/STARsolo/{RECIPE}/Solo.out/GeneFull/raw/UniqueAndMult-EM.mtx",
#         GENEFULL_BC="{OUTDIR}/{SAMPLE}/short_read/STARsolo/{RECIPE}/Solo.out/GeneFull/raw/barcodes.tsv",
#         GENEFULL_FEAT="{OUTDIR}/{SAMPLE}/short_read/STARsolo/{RECIPE}/Solo.out/GeneFull/raw/features.tsv",
#     params:
#         STAR_REF=lambda w: get_STAR_ref(w),
#         STAR_PARAMS=lambda w: get_STAR_extra_params(w),
#         WHITELIST=lambda w: " ".join(get_whitelist(w, return_type="list")),  # space-delimited for multi-barcode
#     resources:
#         mem=megabytes2bytes(config["MEMLIMIT_MB"]),
#         time="2:00:00",
#     threads: config["CORES"]
#     # conda:
#     #     f"{workflow.basedir}/envs/star.yml"
#     priority: 42
#     run:
#         shell(
#             f"""
#             mkdir -p $(dirname {output.BAM})

#             STAR \
#                 --runThreadN {threads} \
#                 --outFileNamePrefix $(dirname {output.BAM})/ \
#                 --outSAMtype BAM SortedByCoordinate \
#                 --limitBAMsortRAM={resources.mem} \
#                 --outSAMattributes NH HI nM AS CR UR CB UB GX GN sS sQ sM \
#                 --readFilesCommand zcat \
#                 --genomeDir {params.STAR_REF} \
#                 --readFilesIn {input.FQS[1]} {input.FQS[0]} \
#                 --outReadsUnmapped Fastx \
#                 --outSAMunmapped Within KeepPairs \
#                 --soloType {params.STAR_PARAMS['STAR_soloType']} {params.STAR_PARAMS['STAR_soloUMI']} {params.STAR_PARAMS['STAR_soloCB']} {params.STAR_PARAMS['STAR_soloAdapter']} {params.STAR_PARAMS['STAR_extra']} \
#                 --soloCBmatchWLtype {params.STAR_PARAMS['STAR_soloCBmatchWLtype']} \
#                 --soloCBwhitelist {params.WHITELIST} \
#                 --soloCellFilter TopCells $(wc -l {input.WHITELIST[0]}) \
#                 --soloUMIfiltering MultiGeneUMI CR \
#                 --soloUMIdedup 1MM_CR \
#                 --soloBarcodeReadLength 0 \
#                 --soloFeatures Gene GeneFull Velocyto \
#                 --soloMultiMappers EM
#             """
#         )


# First pass of STAR alignment
# TODO: add temp()
rule ilmn_3a_STARsolo_firstPass:
    input:
        FQS=lambda w: get_fqs(w, return_type="list", mode="ILMN"),
    output:
        SJ="{OUTDIR}/{SAMPLE}/short_read/STARsolo/{RECIPE}/firstPass/_STARpass1/SJ.out.tab",
    params:
        STAR_REF=lambda w: get_STAR_ref(w),
        STAR_PARAMS=lambda w: get_STAR_extra_params(w),
        STAR_EXTRA=lambda w: get_STAR_extra_params(w)["STAR_extra"],
        WHITELIST=lambda w: " ".join(get_whitelist(w, return_type="list")),  # space-delimited for multi-barcode
    log:
        log="{OUTDIR}/{SAMPLE}/short_read/STARsolo/{RECIPE}/firstPass/pass1.log",
        err="{OUTDIR}/{SAMPLE}/short_read/STARsolo/{RECIPE}/firstPass/pass1.err",
    resources:
        mem=megabytes2bytes(config["MEMLIMIT_MB"]),
        time="2:00:00",
    threads: config["CORES"]
    conda:
        f"{workflow.basedir}/envs/star.yml"
    priority: 42
    shell:
        """
        mkdir -p $(dirname {output.SJ})

        STAR --runThreadN {threads} \
            --outFileNamePrefix $(dirname {output.SJ})/ \
            --genomeDir {params.STAR_REF} \
            --readFilesCommand zcat \
            --readFilesIn {input.FQS[1]} \
            --outSAMtype BAM Unsorted {params.STAR_EXTRA} \
            --twopassMode Basic \
        1> {log.log} \
        2> {log.err}
        """


# Second pass of STAR alignment; includes barcode calling, etc.
# TODO- add sj filtering? (https://www.reneshbedre.com/blog/star-aligner-twopass-mode.html)
rule ilmn_3a_STARsolo_secondPass:
    input:
        FQS=lambda w: get_fqs(w, return_type="list", mode="ILMN"),
        WHITELIST=lambda w: get_whitelist(w, return_type="list"),
        SJ="{OUTDIR}/{SAMPLE}/short_read/STARsolo/{RECIPE}/firstPass/_STARpass1/SJ.out.tab",
    output:
        BAM="{OUTDIR}/{SAMPLE}/short_read/STARsolo/{RECIPE}/Aligned.sortedByCoord.out.bam",
        UNMAPPED1="{OUTDIR}/{SAMPLE}/short_read/STARsolo/{RECIPE}/Unmapped.out.mate1",
        UNMAPPED2="{OUTDIR}/{SAMPLE}/short_read/STARsolo/{RECIPE}/Unmapped.out.mate2",
        VEL=[
            f"{{OUTDIR}}/{{SAMPLE}}/short_read/STARsolo/{{RECIPE}}/Solo.out/Velocyto/raw/{FILE}"
            for FILE in ["spliced.mtx", "barcodes.tsv", "features.tsv"]
        ],
        GENE=[
            f"{{OUTDIR}}/{{SAMPLE}}/short_read/STARsolo/{{RECIPE}}/Solo.out/Gene/raw/{FILE}"
            for FILE in [
                "matrix.mtx",
                "UniqueAndMult-EM.mtx",
                "barcodes.tsv",
                "features.tsv",
            ]
        ],
        GENEFULL=[
            f"{{OUTDIR}}/{{SAMPLE}}/short_read/STARsolo/{{RECIPE}}/Solo.out/GeneFull/raw/{FILE}"
            for FILE in [
                "matrix.mtx",
                "UniqueAndMult-EM.mtx",
                "barcodes.tsv",
                "features.tsv",
            ]
        ],
    params:
        WHITELIST=lambda w: " ".join(get_whitelist(w, return_type="list")),  # space-delimited for multi-barcode
        FQS=lambda w: " ".join(reversed(get_fqs(w, return_type="list", mode="ILMN"))),
        N_CELLS=lambda w: get_n_cells(w),
        STAR_REF=lambda w: get_STAR_ref(w),
        # STAR_PARAMS=lambda w: get_STAR_extra_params(w),
        STAR_PARAMS=lambda w: " ".join(get_STAR_extra_params(w).values()),
        STAR_soloCBmatchWLtype=lambda w: get_recipe_info(w, "STAR_soloCBmatchWLtype"),
        outBAMsortingBinsN=100,  # default: 50; increase to save mem usage
        outBAMsortingThreadN=4,  # default: min(6,–runThreadN); reduce to save mem
    log:
        log="{OUTDIR}/{SAMPLE}/short_read/STARsolo/{RECIPE}/pass2.log",
        err="{OUTDIR}/{SAMPLE}/short_read/STARsolo/{RECIPE}/pass2.err",
    resources:
        # mem=megabytes2bytes(config["MEMLIMIT_MB"]),
        mem=config["MEMLIMIT"],
        time="2:00:00",
    threads: config["CORES"]
    conda:
        f"{workflow.basedir}/envs/star.yml"
    priority: 42
    shell:
        """
        mkdir -p $(dirname {output.BAM})

        STAR --readFilesCommand zcat --readFilesIn {params.FQS} \
            --genomeDir {params.STAR_REF} \
            --runThreadN {threads} \
            --limitBAMsortRAM={resources.mem} \
            --outBAMsortingBinsN {params.outBAMsortingBinsN} \
            --outBAMsortingThreadN {params.outBAMsortingThreadN} \
            --outFileNamePrefix $(dirname {output.BAM})/ \
            --outSAMtype BAM SortedByCoordinate \
            --outSAMattributes NH HI nM AS CR UR CB UB GX GN sS sQ sM \
            --outReadsUnmapped Fastx --outSAMunmapped Within KeepPairs \
            --soloType {params.STAR_PARAMS} \
            --soloCBmatchWLtype {params.STAR_soloCBmatchWLtype} \
            --soloCBwhitelist {params.WHITELIST} \
            --soloCellFilter TopCells {params.N_CELLS} \
            --soloUMIfiltering MultiGeneUMI CR \
            --soloUMIdedup 1MM_CR \
            --soloBarcodeReadLength 0 \
            --soloFeatures Gene GeneFull Velocyto \
            --soloMultiMappers EM \
            --sjdbFileChrStartEnd {input.SJ} \
        1> {log.log} \
        2> {log.err}
        """


#
# --soloType {params.STAR_PARAMS['STAR_soloType']} {params.STAR_PARAMS['STAR_soloUMI']} {params.STAR_PARAMS['STAR_soloCB']} {params.STAR_PARAMS['STAR_soloAdapter']} {params.STAR_PARAMS['STAR_extra']} \


# compress outputs from STAR (count matrices, cell barcodes, and gene lists)
rule ilmn_3a_compress_STAR_outs:
    input:
        VEL_MAT="{OUTDIR}/{SAMPLE}/short_read/STARsolo/{RECIPE}/Solo.out/Velocyto/raw/spliced.mtx",
        VEL_BC="{OUTDIR}/{SAMPLE}/short_read/STARsolo/{RECIPE}/Solo.out/Velocyto/raw/barcodes.tsv",
        VEL_FEAT="{OUTDIR}/{SAMPLE}/short_read/STARsolo/{RECIPE}/Solo.out/Velocyto/raw/features.tsv",
        GENE_MAT="{OUTDIR}/{SAMPLE}/short_read/STARsolo/{RECIPE}/Solo.out/Gene/raw/matrix.mtx",
        GENE_MAT_EM="{OUTDIR}/{SAMPLE}/short_read/STARsolo/{RECIPE}/Solo.out/Gene/raw/UniqueAndMult-EM.mtx",
        GENE_BC="{OUTDIR}/{SAMPLE}/short_read/STARsolo/{RECIPE}/Solo.out/Gene/raw/barcodes.tsv",
        GENE_FEAT="{OUTDIR}/{SAMPLE}/short_read/STARsolo/{RECIPE}/Solo.out/Gene/raw/features.tsv",
        GENEFULL_MAT="{OUTDIR}/{SAMPLE}/short_read/STARsolo/{RECIPE}/Solo.out/GeneFull/raw/matrix.mtx",
        GENEFULL_MAT_EM="{OUTDIR}/{SAMPLE}/short_read/STARsolo/{RECIPE}/Solo.out/GeneFull/raw/UniqueAndMult-EM.mtx",
        GENEFULL_BC="{OUTDIR}/{SAMPLE}/short_read/STARsolo/{RECIPE}/Solo.out/GeneFull/raw/barcodes.tsv",
        GENEFULL_FEAT="{OUTDIR}/{SAMPLE}/short_read/STARsolo/{RECIPE}/Solo.out/GeneFull/raw/features.tsv",
    output:
        VEL_MAT="{OUTDIR}/{SAMPLE}/short_read/STARsolo/{RECIPE}/Solo.out/Velocyto/raw/spliced.mtx.gz",
        VEL_BC="{OUTDIR}/{SAMPLE}/short_read/STARsolo/{RECIPE}/Solo.out/Velocyto/raw/barcodes.tsv.gz",
        VEL_FEAT="{OUTDIR}/{SAMPLE}/short_read/STARsolo/{RECIPE}/Solo.out/Velocyto/raw/features.tsv.gz",
        GENE_MAT="{OUTDIR}/{SAMPLE}/short_read/STARsolo/{RECIPE}/Solo.out/Gene/raw/matrix.mtx.gz",
        GENE_MAT_EM="{OUTDIR}/{SAMPLE}/short_read/STARsolo/{RECIPE}/Solo.out/Gene/raw/UniqueAndMult-EM.mtx.gz",
        GENE_BC="{OUTDIR}/{SAMPLE}/short_read/STARsolo/{RECIPE}/Solo.out/Gene/raw/barcodes.tsv.gz",
        GENE_FEAT="{OUTDIR}/{SAMPLE}/short_read/STARsolo/{RECIPE}/Solo.out/Gene/raw/features.tsv.gz",
        GENEFULL_MAT="{OUTDIR}/{SAMPLE}/short_read/STARsolo/{RECIPE}/Solo.out/GeneFull/raw/matrix.mtx.gz",
        GENEFULL_MAT_EM="{OUTDIR}/{SAMPLE}/short_read/STARsolo/{RECIPE}/Solo.out/GeneFull/raw/UniqueAndMult-EM.mtx.gz",
        GENEFULL_BC="{OUTDIR}/{SAMPLE}/short_read/STARsolo/{RECIPE}/Solo.out/GeneFull/raw/barcodes.tsv.gz",
        GENEFULL_FEAT="{OUTDIR}/{SAMPLE}/short_read/STARsolo/{RECIPE}/Solo.out/GeneFull/raw/features.tsv.gz",
    params:
        VELDIR=directory(
            "{OUTDIR}/{SAMPLE}/short_read/STARsolo/{RECIPE}/Solo.out/Velocyto"
        ),
        GENEDIR=directory(
            "{OUTDIR}/{SAMPLE}/short_read/STARsolo/{RECIPE}/Solo.out/Gene"
        ),
        GENEFULLDIR=directory(
            "{OUTDIR}/{SAMPLE}/short_read/STARsolo/{RECIPE}/Solo.out/GeneFull"
        ),
    # resources:
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
            pigz -p{threads} -f \
                {OUTDIR}/{wildcards.SAMPLE}/short_read/STARsolo/{wildcards.RECIPE}/*/*/*/*.tsv \
                {OUTDIR}/{wildcards.SAMPLE}/short_read/STARsolo/{wildcards.RECIPE}/*/*/*/*.mtx
            """
        )
