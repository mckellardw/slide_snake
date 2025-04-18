#############################################
## Unmapped read analyses
#############################################


# Run fastqc on unmapped reads; switch names because of STAR weirdness
rule ilmn_3b_fastqc_unmapped:
    input:
        UNMAPPED1="{OUTDIR}/{SAMPLE}/short_read/STARsolo/{RECIPE}/Unmapped.out.mate1",
        UNMAPPED2="{OUTDIR}/{SAMPLE}/short_read/STARsolo/{RECIPE}/Unmapped.out.mate2",
    output:
        UNMAPPED1_FQ="{OUTDIR}/{SAMPLE}/short_read/STARsolo/{RECIPE}/Unmapped.out.mate1.fq.gz",
        UNMAPPED2_FQ="{OUTDIR}/{SAMPLE}/short_read/STARsolo/{RECIPE}/Unmapped.out.mate2.fq.gz",
        fastqcDir=directory("{OUTDIR}/{SAMPLE}/short_read/fastqc/unmapped/{RECIPE}"),
    params:
        FASTQC_ADAPTERS=config["FASTQC_ADAPTERS"],
    log:
        log="{OUTDIR}/{SAMPLE}/short_read/fastqc/unmapped/{RECIPE}/fastqc_unmapped.log",
    threads: config["CORES"]
    conda:
        f"{workflow.basedir}/envs/fastqc.yml"
    shell:
        """        
        mv {input.UNMAPPED1} {input.UNMAPPED2}.fq
        mv {input.UNMAPPED2} {input.UNMAPPED1}.fq

        pigz -p{threads} -f {input.UNMAPPED1}.fq {input.UNMAPPED2}.fq

        mkdir -p {output.fastqcDir}

        fastqc \
            -o {output.fastqcDir} \
            -t {threads} \
            -a {params.FASTQC_ADAPTERS} \
            {output.UNMAPPED1_FQ} {output.UNMAPPED2_FQ} \
            > {log.log} 2>&1
        """


# TODO - de novo assembly of unmapped reads?
# TODO - microbiome analysis of unmapped reads?
