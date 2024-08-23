#TODO
# github: https://github.com/COMBINE-lab/oarfish
# documentation: 

rule ont_kallisto_lr:
    input:
        FQ=lambda w: get_fqs(w, return_type="list", mode="ONT"),
        BC="{OUTDIR}/{SAMPLE}/bc/whitelist.txt",
    output:
        SAM_TMP="{OUTDIR}/{SAMPLE}/ont/kb_lr/{RECIPE}/",
    params:
        EXTRA_FLAGS=lambda wildcards: RECIPE_SHEET["minimap2.extra"][wildcards.RECIPE],
        ref=config["REF_GENOME_FASTA"],
        chrom_sizes=config["REF_CHROM_SIZES"],
        bed=config["REF_GENES_BED"],
        # flags=config["RESOURCES_MM2_FLAGS"],
    log:
        log="{OUTDIR}/{SAMPLE}/ont/kb_lr/{RECIPE}/kallisto_lr.log",
    resources:
        mem="128G",
    threads: config["CORES"]
    conda:
        f"{workflow.basedir}/envs/kb_lr.yml"
    shell:
        f"""
        bash scripts/bash/kb_lr.sh \
            --outdir $(dirname {output.BUS}) \
            --kb_idx {KB_IDX} \
            --whitelist {input.BC} \
            --chemistry {KB_X} \
            --log {log.log} \
            --threads {resources.threads} \
            --memlimit {params.mem} \
            --r1fq {input.FQS[0]} \
            --r2fq {input.FQS[1]}
        """



rule ont_bus2mat_lr:
    input:
        BUS="{OUTDIR}/{SAMPLE}/kb/{RECIPE}/output.corrected.bus",
        TRANSCRIPTS="{OUTDIR}/{SAMPLE}/kb/{RECIPE}/transcripts.txt",
        ECMAP="{OUTDIR}/{SAMPLE}/kb/{RECIPE}/matrix.ec",
    output:
        BCS="{OUTDIR}/{SAMPLE}/kb/{RECIPE}/raw/output.barcodes.txt",
        GENES="{OUTDIR}/{SAMPLE}/kb/{RECIPE}/raw/output.genes.txt",
        MAT="{OUTDIR}/{SAMPLE}/kb/{RECIPE}/raw/output.mtx",
    log:
        log="{OUTDIR}/{SAMPLE}/kb/{RECIPE}/raw/bustools_count.log",
    threads: 1
    run:
        KB_T2G = T2G_DICT[wildcards.SAMPLE]

        shell(
            f"""
            mkdir -p $(dirname {output.MAT})

            bustools count \
                --output $(dirname {output.MAT})/ \
                --genemap {KB_T2G} \
                --ecmap {input.ECMAP} \
                --txnames {input.TRANSCRIPTS} \
                --genecounts \
                --umi-gene \
                --em \
                {input.BUS} \
            2> {log.log}
            """
        )


# gzip the count matrix, etc.
rule compress_kb_outs:
    input:
        BCS="{OUTDIR}/{SAMPLE}/kb/{RECIPE}/raw/output.barcodes.txt",
        GENES="{OUTDIR}/{SAMPLE}/kb/{RECIPE}/raw/output.genes.txt",
        MAT="{OUTDIR}/{SAMPLE}/kb/{RECIPE}/raw/output.mtx",
    output:
        BCS="{OUTDIR}/{SAMPLE}/kb/{RECIPE}/raw/output.barcodes.txt.gz",
        GENES="{OUTDIR}/{SAMPLE}/kb/{RECIPE}/raw/output.genes.txt.gz",
        MAT="{OUTDIR}/{SAMPLE}/kb/{RECIPE}/raw/output.mtx.gz",
    threads: config["CORES"]
    shell:
        f"""
        pigz -p{resources.threads} {input.BCS} {input.GENES} {input.MAT}
        """
