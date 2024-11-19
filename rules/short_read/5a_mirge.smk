#############################################
## miRge3.0 analysis
#############################################
#TODO- add rule to filter out longer reads for faster smRNA analysis?


#TODO- remove this when miRge3 updates to allow '.fq.gz' as inputs...
rule copy_R2_fq_for_mirge:
    input:
        FQ=lambda w: get_fqs(w, return_type="list", mode="ILMN")[1],
    output:
        FQ=temp("{OUTDIR}/{SAMPLE}/short_read/miRge_bulk/{RECIPE}/tmp/R2.fastq.gz"),
    shell:
        """
        cp {input.FQ} {output.FQ}
        """


# Source: https://mirge3.readthedocs.io/en/latest/quick_start.html
## Note- `--outDirNam` is a hidden argument for miRge3 that allows direct naming of the output directory
# TODO update this code...
rule miRge3_pseudobulk:
    input:
        FQ="{OUTDIR}/{SAMPLE}/short_read/miRge_bulk/{RECIPE}/tmp/R2.fastq.gz",
    output:
        MIRGE_HTML="{OUTDIR}/{SAMPLE}/short_read/miRge_bulk/{RECIPE}/annotation.report.html",
    params:
        MIRGE_LIB=os.path.abspath(config["MIRGE_LIB"]),
        SPECIES=lambda wildcards: SAMPLE_SHEET["species"][wildcards.SAMPLE],
        MIN_LENGTH = 12,
    threads: config["CORES"]
    resources:
        mem="64G",
    conda:
        f"{workflow.basedir}/envs/mirge3.yml"
    shell:
        """
        mkdir -p $(dirname {output.MIRGE_HTML})
        cd $(dirname {output.MIRGE_HTML})

        miRge3.0 \
            -s {input.FQ} \
            -lib {params.MIRGE_LIB} \
            -on {params.SPECIES} \
            -db mirbase \
            --minimum-length {params.MIN_LENGTH} \
            --outDirName ./ \
            --threads {threads} \
            --minReadCounts 1 \
            --gff-out \
            --novel-miRNA \
            --AtoI
        """
 
        #TODO- extra flags for human:
            # --miREC {EXTRA_FLAGS}
        #     # human-only settings
        # if params.SPECIES == "human":
        #     EXTRA_FLAGS = "--tRNA-frag"
        # else:
        #     EXTRA_FLAGS = ""
