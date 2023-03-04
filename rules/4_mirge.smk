#############################################
## miRge3.0 analysis
#############################################

#TODO- add rule to filter out longer reads for faster smRNA analysis?


#Source: https://mirge3.readthedocs.io/en/latest/quick_start.html
## Note- `--outDirNam` is a hidden argument for miRge3 that allows direct naming of the output directory
rule miRge3_pseudobulk:
    input:
        FINAL_R2_FQ = '{OUTDIR}/{sample}/tmp/{sample}_R2_final.fq.gz'
    output:
        MIRGE_DIR = directory('{OUTDIR}/{sample}/miRge_bulk'),
        MIRGE_HTML = '{OUTDIR}/{sample}/miRge_bulk/annotation.report.html'
        # MIRGE_CHECK = '{OUTDIR}/{sample}/miRge_check.txt'
    params:
        MIRGE_LIB = config['MIRGE_LIB'],
        # SPECIES = config['SPECIES'],
        # UMIlen = config['UMIlen'],
        MEMLIMIT = config['MEMLIMIT']
    threads:
        config['CORES']
    run:
        SPECIES = SPECIES_DICT[wildcards.sample] #TODO

        shell(
            f"""
            zcat {input.FINAL_R2_FQ} > {OUTDIR}/{wildcards.sample}/tmp/tmp_R2.fastq

            {MIRGE_EXEC} \
            -s {OUTDIR}/{wildcards.sample}/tmp/tmp_R2.fastq \
            -lib {params.MIRGE_LIB} \
            -on {SPECIES} \
            -db mirbase \
            --nextseq-trim 1 \
            --minimum-length 12 \
            --outDirName {output.MIRGE_DIR} \
            --threads {threads} \
            -gff -nmir -ai             

            rm {OUTDIR}/{wildcards.sample}/tmp/tmp_R2.fastq
            """
            # mkdir -p {output.MIRGE_DIR}
            # mv --backup=numbered {OUTDIR}/{wildcards.sample}/miRge.*/* {output.MIRGE_DIR}
            # rmdir {OUTDIR}/{wildcards.sample}/miRge.*
            
        )
        # && touch {output.MIRGE_CHECK}
# 
# -a illumina \
# -trf

# TODO- miRge across cells/spots/barcodes