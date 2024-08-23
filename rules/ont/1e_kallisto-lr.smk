#TODO
# github: https://github.com/COMBINE-lab/oarfish
# documentation: 

rule ont_kallisto_lr:
    input:
        FQ=lambda w: get_fqs(w, return_type="list", mode="ONT"),
        BC="{OUTDIR}/{SAMPLE}/bc/whitelist.txt",
    output:
        BUS=temp("{OUTDIR}/{SAMPLE}/ont/kb_lr/{RECIPE}/output.bus"),
        BUS_CORRECTED=temp("{OUTDIR}/{SAMPLE}/ont/kb_lr/{RECIPE}/output.corrected.bus"),
        TRANSCRIPTS="{OUTDIR}/{SAMPLE}/ont/kb_lr/{RECIPE}/transcripts.txt",
        ECMAP=temp("{OUTDIR}/{SAMPLE}/ont/kb_lr/{RECIPE}/matrix.ec"),
    params:
        KB_IDX=lambda wildcards: IDX_DICT[wildcards.SAMPLE],
        KB_X=lambda wildcards: RECIPE_SHEET["kb.x"][wildcards.RECIPE],
        KB_EXTRA=lambda wildcards: RECIPE_SHEET["kb.extra"][wildcards.RECIPE],
        KB_T2G=lambda wildcards: T2G_DICT[wildcards.SAMPLE],
        N_READS_SUMMARY=1000000,  # number of reads to use for summary stats
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
            --kb_idx {params.KB_IDX} \
            --whitelist {input.BC} \
            --chemistry {params.KB_X} \
            --log {log.log} \
            --threads {threads} \
            --memlimit {resources.mem} \
            --r1fq {input.FQS[0]} \
            --r2fq {input.FQS[1]}
        """



rule ont_bus2mat_lr:
    input:
        BUS="{OUTDIR}/{SAMPLE}/ont/kb_lr/{RECIPE}/output.corrected.bus",
        TRANSCRIPTS="{OUTDIR}/{SAMPLE}/ont/kb_lr/{RECIPE}/transcripts.txt",
        ECMAP="{OUTDIR}/{SAMPLE}/ont/kb_lr/{RECIPE}/matrix.ec",
    output:
        BCS="{OUTDIR}/{SAMPLE}/ont/kb_lr{RECIPE}/raw/output.barcodes.txt",
        GENES="{OUTDIR}/{SAMPLE}/ont/kb_lr{RECIPE}/raw/output.genes.txt",
        MAT="{OUTDIR}/{SAMPLE}/ont/kb_lr{RECIPE}/raw/output.mtx",
    log:
        log="{OUTDIR}/{SAMPLE}/ont/kb_lr{RECIPE}/raw/bustools_count.log",
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
        BCS="{OUTDIR}/{SAMPLE}/ont/kb_lr{RECIPE}/raw/output.barcodes.txt",
        GENES="{OUTDIR}/{SAMPLE}/ont/kb_lr{RECIPE}/raw/output.genes.txt",
        MAT="{OUTDIR}/{SAMPLE}/ont/kb_lr{RECIPE}/raw/output.mtx",
    output:
        BCS="{OUTDIR}/{SAMPLE}/ont/kb_lr{RECIPE}/raw/output.barcodes.txt.gz",
        GENES="{OUTDIR}/{SAMPLE}/ont/kb_lr{RECIPE}/raw/output.genes.txt.gz",
        MAT="{OUTDIR}/{SAMPLE}/ont/kb_lr{RECIPE}/raw/output.mtx.gz",
    threads: config["CORES"]
    shell:
        f"""
        pigz -p{resources.threads} {input.BCS} {input.GENES} {input.MAT}
        """
