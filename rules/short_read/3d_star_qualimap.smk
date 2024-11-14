# ALignment qc with qualimap
## link: https://qualimap.conesalab.org/


## qualimap on aligned reads (not deduplicated)
rule ilmn_3d_qualimapQC_STAR:
    input:
        BAM="{OUTDIR}/{SAMPLE}/short_read/STARsolo/{RECIPE}/Aligned.sortedByCoord.out.bam",
        BAI="{OUTDIR}/{SAMPLE}/short_read/STARsolo/{RECIPE}/Aligned.sortedByCoord.out.bam.bai",
    output:
        TXT="{OUTDIR}/{SAMPLE}/short_read/qualimap/STAR/{RECIPE}/raw/rnaseq_qc_results.txt",
        HTML="{OUTDIR}/{SAMPLE}/short_read/qualimap/STAR/{RECIPE}/raw/qualimapReport.html",
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

## qualimap on deduplicated/aligned reads
rule ilmn_3d_qualimapQC_dedup_STAR:
    input:
        BAM="{OUTDIR}/{SAMPLE}/short_read/STARsolo/{RECIPE}/Aligned.sortedByCoord.out.dedup.bam",
        BAI="{OUTDIR}/{SAMPLE}/short_read/STARsolo/{RECIPE}/Aligned.sortedByCoord.out.dedup.bam.bai",
    output:
        TXT="{OUTDIR}/{SAMPLE}/short_read/qualimap/STAR/{RECIPE}/dedup/rnaseq_qc_results.txt",
        HTML="{OUTDIR}/{SAMPLE}/short_read/qualimap/STAR/{RECIPE}/dedup/qualimapReport.html",
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


rule ilmn_3d_qualimap_summary2csv_STAR:
    input:
        TXT="{OUTDIR}/{SAMPLE}/short_read/qualimap/STAR/{RECIPE}/{DEDUP}/rnaseq_qc_results.txt",
    output:
        CSV="{OUTDIR}/{SAMPLE}/short_read/qualimap/STAR/{RECIPE}/{DEDUP}/rnaseq_qc_result.csv",
    threads: 1
    shell:
        """
        python scripts/py/qualimap_summary2csv.py {input.TXT} {output.CSV}
        """
