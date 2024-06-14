
rule ont_qualimap_minimap2:
    input:
        BAM="{OUTDIR}/{SAMPLE}/ont/minimap2/{RECIPE}/sorted_gn_cb_ub.bam",
        BAI="{OUTDIR}/{SAMPLE}/ont/minimap2/{RECIPE}/sorted_gn_cb_ub.bam.bai",
    output:
        TXT="{OUTDIR}/{SAMPLE}/qualimap/ont/minimap2/{RECIPE}/rnaseq_qc_results.txt",
        HTML="{OUTDIR}/{SAMPLE}/qualimap/ont/minimap2/{RECIPE}/qualimapReport.html",
    params:
        GENES_GTF=lambda wildcards: GTF_DICT[wildcards.SAMPLE],
    threads: 1
    resources:
        mem="32G",
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
            -outformat html
        """


rule ont_qualimap_STAR:
    input:
        BAM="{OUTDIR}/{SAMPLE}/STARsolo/ont/{RECIPE}/Aligned.sortedByCoord.out.bam",
        BAI="{OUTDIR}/{SAMPLE}/STARsolo/ont/{RECIPE}/Aligned.sortedByCoord.out.bam.bai",
    output:
        TXT="{OUTDIR}/{SAMPLE}/qualimap/ont/STARsolo/{RECIPE}/rnaseq_qc_results.txt",
        HTML="{OUTDIR}/{SAMPLE}/qualimap/ont/STARsolo/{RECIPE}/qualimapReport.html",
    params:
        GENES_GTF = lambda wildcards: GTF_DICT[wildcards.SAMPLE]
    threads: 1
    resources:
        mem="32G",
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
            -outformat html
        """


rule ont_qualimap_summary2csv:
    input:
        TXT="{OUTDIR}/{SAMPLE}/qualimap/ont/{TOOL}/rnaseq_qc_results.txt",
    output:
        CSV="{OUTDIR}/{SAMPLE}/qualimap/ont/{TOOL}/rnaseq_qc_results.csv",
    threads: 1
    run:
        shell(
            f"""
            python scripts/py/qualimap_summary2csv.py {input.TXT} {output.CSV}
            """
        )
