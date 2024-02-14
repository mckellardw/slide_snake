#############################################
## miRge3.0 analysis
#############################################

#TODO- add rule to filter out longer reads for faster smRNA analysis?

#TODO- remove this when miRge3 updates to allow '.fq.gz' as inputs...
rule copy_fq_for_mirge:
    input:
        R2_FQ = '{OUTDIR}/{SAMPLE}/tmp/cut_R2.fq.gz',
        R2_FQ_TWICE_CUT = '{OUTDIR}/{SAMPLE}/tmp/twiceCut_R2.fq.gz',
        R2_FQ_STAR_FILTERED = '{OUTDIR}/{SAMPLE}/rRNA/STARsolo/final_filtered_R2.fq.gz',
        R2_FQ_BWA_FILTERED  = '{OUTDIR}/{SAMPLE}/rRNA/bwa/final_filtered_R2.fq.gz',
    output:
        R2_FQ = temp('{OUTDIR}/{SAMPLE}/tmp/cut_R2.fastq.gz'),
        R2_FQ_TWICE_CUT = temp('{OUTDIR}/{SAMPLE}/tmp/twiceCut_R2.fastq.gz'),
        R2_FQ_STAR_FILTERED = temp('{OUTDIR}/{SAMPLE}/rRNA/STARsolo/final_filtered_R2.fastq.gz'),
        R2_FQ_BWA_FILTERED  = temp('{OUTDIR}/{SAMPLE}/rRNA/bwa/final_filtered_R2.fastq.gz'),
    run:
        shell(
            f"""
            cp {input.R2_FQ} {output.R2_FQ}
            cp {input.R2_FQ_TWICE_CUT} {output.R2_FQ_TWICE_CUT}
            cp {input.R2_FQ_STAR_FILTERED} {output.R2_FQ_STAR_FILTERED}
            cp {input.R2_FQ_BWA_FILTERED} {output.R2_FQ_BWA_FILTERED}
            """
        )


#Source: https://mirge3.readthedocs.io/en/latest/quick_start.html
## Note- `--outDirNam` is a hidden argument for miRge3 that allows direct naming of the output directory
rule miRge3_pseudobulk:
    input:
        R2_FQ = '{OUTDIR}/{SAMPLE}/tmp/cut_R2.fastq.gz',
        R2_FQ_TWICE_CUT = '{OUTDIR}/{SAMPLE}/tmp/twiceCut_R2.fastq.gz',
        R2_FQ_STAR_FILTERED = '{OUTDIR}/{SAMPLE}/rRNA/STARsolo/final_filtered_R2.fastq.gz',
        R2_FQ_BWA_FILTERED  = '{OUTDIR}/{SAMPLE}/rRNA/bwa/final_filtered_R2.fastq.gz',
    output:
        # GUNZIP_R2_FQ = temp('{OUTDIR}/{SAMPLE}/tmp/twiceCut_R2.fq'),
        MIRGE_DIR = directory('{OUTDIR}/{SAMPLE}/miRge_bulk/{RECIPE}'),
        MIRGE_HTML = '{OUTDIR}/{SAMPLE}/miRge_bulk/{RECIPE}/annotation.report.html'
    params:
        MIRGE_LIB = config['MIRGE_LIB'],
        # SPECIES = config['SPECIES'],
        # UMIlen = config['UMIlen'],
        MEMLIMIT = config['MEMLIMIT']
    threads:
        config['CORES']
    run:
        from os import path

        SPECIES = SPECIES_DICT[wildcards.SAMPLE]

        # recipe = RECIPE_DICT[wildcards.SAMPLE]
        recipe = wildcards.RECIPE
        
        # Select input reads based on alignment recipe
        if "rRNA.STAR" in recipe: # Use trimmed & STAR-rRNA-filtered .fq's
            R2 = input.R2_FQ_STAR_FILTERED
        elif "rRNA.bwa" in recipe: #TODO Use trimmed & bwa-rRNA-filtered .fq's
            R2 = input.R2_FQ_BWA_FILTERED
        elif "rRNA" not in recipe: # just trimmed .fq's
            # R2 = input.R2_FQ
            R2 = input.R2_FQ_TWICE_CUT
        else:
            print("I just don't know what to do with myself...")

        # human-only settings
        if SPECIES == "human":
            EXTRA_FLAGS = "--tRNA-frag"
        else:
            EXTRA_FLAGS = ""

        MIRGE_LIB_ABS = path.abspath(params.MIRGE_LIB)

        # zcat {R2} > {R2.strip('.gz')}
        shell(
            f"""
            mkdir -p {output.MIRGE_DIR}
            cd {output.MIRGE_DIR}

            {EXEC['MIRGE']} \
                -s {R2} \
                -lib {MIRGE_LIB_ABS} \
                -on {SPECIES} \
                -db mirbase \
                --minimum-length 12 \
                --outDirName ./ \
                --threads {threads} \
                --minReadCounts 1 \
                --gff-out \
                --novel-miRNA \
                --AtoI \
                --miREC {EXTRA_FLAGS}
            """            
        )
# -a illumina \
# -trf

# TODO split bame across barcodes
# TODO convert split bams to fqs
# TODO- miRge across cells/spots/barcodes