#############################################
## miRge3.0 analysis
#############################################

#TODO- add rule to filter out longer reads for faster smRNA analysis?


#Source: https://mirge3.readthedocs.io/en/latest/quick_start.html
## Note- `--outDirNam` is a hidden argument for miRge3 that allows direct naming of the output directory
rule miRge3_pseudobulk:
    input:
        R2_FQ = '{OUTDIR}/{SAMPLE}/tmp/twiceCut_R2.fq.gz'
    output:
        GIUNZIP_R2_FQ = temp('{OUTDIR}/{SAMPLE}/tmp/twiceCut_R2.fq'),
        MIRGE_DIR = directory('{OUTDIR}/{SAMPLE}/miRge_bulk'),
        MIRGE_HTML = '{OUTDIR}/{SAMPLE}/miRge_bulk/annotation.report.html'
        # MIRGE_CHECK = '{OUTDIR}/{SAMPLE}/miRge_check.txt'
    params:
        MIRGE_LIB = config['MIRGE_LIB'],
        # SPECIES = config['SPECIES'],
        # UMIlen = config['UMIlen'],
        MEMLIMIT = config['MEMLIMIT']
    threads:
        config['CORES']
    run:
        SPECIES = SPECIES_DICT[wildcards.SAMPLE] #TODO

        shell(
            f"""
            zcat {input.R2_FQ} > {input.R2_FQ.strip('.gz')}

            {EXEC['MIRGE']} \
                -s {input.R2_FQ.strip('.gz')} \
                -lib {params.MIRGE_LIB} \
                -on {SPECIES} \
                -db mirbase \
                --nextseq-trim 1 \
                --minimum-length 12 \
                --outDirName {output.MIRGE_DIR} \
                --threads {threads} \
                -gff -nmir -ai             
            """            
        )
# -a illumina \
# -trf

# TODO split bame across barcodes
# TODO convert split bams to fqs
# TODO- miRge across cells/spots/barcodes