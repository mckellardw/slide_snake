# rRNA removeal w/ RiboDetector
## Source: https://github.com/hzi-bifo/RiboDetector
## Paper: https://doi.org/10.1093/nar/gkac112

rule ribodetector:
    input:
        R1_FQ="{OUTDIR}/{SAMPLE}/tmp/twiceCut_R1.fq.gz",
        R2_FQ="{OUTDIR}/{SAMPLE}/tmp/twiceCut_R2.fq.gz",
    output:
        R1_FQ="{OUTDIR}/{SAMPLE}/rRNA/ribodetector/noRibo_R1.fq",
        R2_FQ="{OUTDIR}/{SAMPLE}/rRNA/ribodetector/noRibo_R2.fq",
        RIBO_R1_FQ="{OUTDIR}/{SAMPLE}/rRNA/ribodetector/noRibo_R1.fq",
        RIBO_R2_FQ="{OUTDIR}/{SAMPLE}/rRNA/ribodetector/noRibo_R2.fq",
    params:
        CHUNK_SIZE=256,
        MIN_ALIGNSCORE=40,
    log:
        log="{OUTDIR}/{SAMPLE}/rRNA/ribodetector/ribodetector.log",
    threads: config["CORES"]
    conda:
        f"{workflow.basedir}/envs/ribodetector.yml"
    shell:
        """
        seqkit stats {input.R2_FQ} > {log.log}
        READ_LEN=$(cat {log.log} | tail -n 1 | awk '{{print $7}}')

        ribodetector_cpu \
            --threads {threads} \
            --len ${READ_LEN} \
            --input {input.R1_FQ} {input.R2_FQ} \
            --ensure rrna \
            --chunk_size {params.CHUNK_SIZE} \
            --rrna {output.R1_FQ} {output.R2_FQ} \
            --output {output.R1_FQ} {output.R2_FQ} \
            --log {log.log}
        """

rule compress_ribodetector:
    input:
        FQ="{OUTDIR}/{SAMPLE}/rRNA/ribodetector/{READ}.fq",
    output:
        FQ="{OUTDIR}/{SAMPLE}/rRNA/ribodetector/{READ}.fq.gz",
    threads: config["CORES"]
    run:
        shell(
            f"""
            {EXEC['PIGZ']} -p{threads} {input.FQ}
            """
        )