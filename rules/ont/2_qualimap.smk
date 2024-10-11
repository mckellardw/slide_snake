
rule ont_2_qualimap_minimap2:
    input:
        BAM="{OUTDIR}/{SAMPLE}/ont/minimap2/{RECIPE}/sorted_gn_cb_ub.bam",
        BAI="{OUTDIR}/{SAMPLE}/ont/minimap2/{RECIPE}/sorted_gn_cb_ub.bam.bai",
    output:
        TXT="{OUTDIR}/{SAMPLE}/ont/qualimap/minimap2/{RECIPE}/rnaseq_qc_results.txt",
        HTML="{OUTDIR}/{SAMPLE}/ont/qualimap/minimap2/{RECIPE}/qualimapReport.html",
    params:
        GENES_GTF=lambda wildcards: SAMPLE_SHEET["genes_gtf"][wildcards.SAMPLE],
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
            -outformat html
        """


rule ont_2_qualimap_STAR:
    input:
        BAM="{OUTDIR}/{SAMPLE}/STARsolo/ont/{RECIPE}/Aligned.sortedByCoord.out.bam",
        BAI="{OUTDIR}/{SAMPLE}/STARsolo/ont/{RECIPE}/Aligned.sortedByCoord.out.bam.bai",
    output:
        TXT="{OUTDIR}/{SAMPLE}/ont/qualimap/STARsolo/{RECIPE}/rnaseq_qc_results.txt",
        HTML="{OUTDIR}/{SAMPLE}/ont/qualimap/STARsolo/{RECIPE}/qualimapReport.html",
    params:
        GENES_GTF=lambda wildcards: SAMPLE_SHEET["genes_gtf"][wildcards.SAMPLE],
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
            -outformat html
        """


rule ont_2_qualimap_summary2csv:
    input:
        TXT="{OUTDIR}/{SAMPLE}/ont/qualimap/{TOOL}/rnaseq_qc_results.txt",
    output:
        CSV="{OUTDIR}/{SAMPLE}/ont/qualimap/{TOOL}/rnaseq_qc_results.csv",
    resources:
        mem="8G",
    threads: 1
    shell:
        """
        python scripts/py/qualimap_summary2csv.py {input.TXT} {output.CSV}
        """
