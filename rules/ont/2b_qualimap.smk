# ALignment qc with qualimap
## link: https://qualimap.conesalab.org/


# Qualimap QC on alignment outputs
rule ont_2b_qualimap:
    input:
        BAM="{OUTDIR}/{SAMPLE}/ont/{REF}/{RECIPE}/sorted_filtered_gn_cb_ub.bam",
        BAI="{OUTDIR}/{SAMPLE}/ont/{REF}/{RECIPE}/sorted_filtered_gn_cb_ub.bam.bai",
    output:
        TXT="{OUTDIR}/{SAMPLE}/ont/qualimap/{REF}/{RECIPE}/rnaseq/rnaseq_qc_results.txt",
        PDF="{OUTDIR}/{SAMPLE}/ont/qualimap/{REF}/{RECIPE}/rnaseq/report.pdf",
    params:
        GENES_GTF=lambda wildcards: SAMPLE_SHEET["genes_gtf"][wildcards.SAMPLE],
    log:
        log="{OUTDIR}/{SAMPLE}/ont/qualimap/{REF}/{RECIPE}/rnaseq/rnaseq.log",
        err="{OUTDIR}/{SAMPLE}/ont/qualimap/{REF}/{RECIPE}/rnaseq/rnaseq.err",
    resources:
        mem="32G",
    threads: 1
    conda:
        f"{workflow.basedir}/envs/qualimap.yml"
    shell:
        """
        mkdir -p $(dirname {output.TXT})

        qualimap rnaseq \
            -bam {input.BAM} \
            -gtf {params.GENES_GTF} \
            --sequencing-protocol strand-specific-forward \
            --sorted \
            --java-mem-size={resources.mem} \
            -outdir $(dirname {output.TXT}) \
            -outformat pdf \
        1> {log.log} \
        2> {log.err}
        """


# Convert the unfortunately formatted qc results from qualimap into a readable format
rule ont_2b_qualimap_readqc_summary2csv:
    input:
        TXT="{OUTDIR}/{SAMPLE}/ont/qualimap/{TOOL}/rnaseq_qc_results.txt",
    output:
        CSV="{OUTDIR}/{SAMPLE}/ont/qualimap/{TOOL}/rnaseq_qc_results.csv",
    resources:
        mem="8G",
    threads: 1
    shell:
        """
        python scripts/py/qualimap_summary2csv.py \
            --input {input.TXT} \
            --output {output.CSV}
        """


# Qualimap BAM QC on alignment outputs
rule ont_2b_qualimap_bamqc:
    input:
        BAM="{OUTDIR}/{SAMPLE}/ont/{REF}/{RECIPE}/sorted_filtered_gn_cb_ub.bam",
        BAI="{OUTDIR}/{SAMPLE}/ont/{REF}/{RECIPE}/sorted_filtered_gn_cb_ub.bam.bai",
    output:
        REPORT="{OUTDIR}/{SAMPLE}/ont/qualimap/{REF}/{RECIPE}/bamqc/report.pdf",
        CSV="{OUTDIR}/{SAMPLE}/ont/qualimap/{REF}/{RECIPE}/bamqc/genome_results.txt",
    log:
        log="{OUTDIR}/{SAMPLE}/ont/qualimap/{REF}/{RECIPE}/bamqc/bamqc.log",
        err="{OUTDIR}/{SAMPLE}/ont/qualimap/{REF}/{RECIPE}/bamqc/bamqc.err",
    resources:
        mem="16G",
    threads: 1
    conda:
        f"{workflow.basedir}/envs/qualimap.yml"
    shell:
        """
        mkdir -p $(dirname {output.REPORT})

        qualimap bamqc \
            -bam {input.BAM} \
            --java-mem-size={resources.mem} \
            -outdir $(dirname {output.REPORT}) \
            -outformat pdf \
        1> {log.log} \
        2> {log.err}
        """


# Convert the unfortunately formatted qc results from qualimap into a readable format
rule ont_2b_qualimap_bamqc_summary2csv:
    input:
        TXT="{OUTDIR}/{SAMPLE}/ont/qualimap/{TOOL}/genome_results.txt",
    output:
        CSV="{OUTDIR}/{SAMPLE}/ont/qualimap/{TOOL}/genome_results.csv",
    resources:
        mem="8G",
    threads: 1
    shell:
        """
        python scripts/py/qualimap_summary2csv.py \
            --input {input.TXT} \
            --output {output.CSV}
        """
