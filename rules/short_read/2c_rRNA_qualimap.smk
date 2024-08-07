# Qualimap on bwa output
## qualimap on deduplicated/aligned reads
rule qualimapQC_rRNA_bwa:
    input:
        BAM="{OUTDIR}/{SAMPLE}/rRNA/bwa/aligned_sorted.bam",
        BAI="{OUTDIR}/{SAMPLE}/rRNA/bwa/aligned_sorted.bam.bai",
    output:
        TXT="{OUTDIR}/{SAMPLE}/qualimap/rRNA/bwa/rnaseq_qc_results.txt",
        HTML="{OUTDIR}/{SAMPLE}/qualimap/rRNA/bwa/qualimapReport.html",
    params:
        GENES_GTF=lambda wildcards: SAMPLE_SHEET["rRNA_gtf"][wildcards.SAMPLE],
    log:
        log="{OUTDIR}/{SAMPLE}/qualimap/rRNA/bwa/rnaseq_qc.log",
    resources:
        mem="32G",
        threads=1,
    conda:
        f"{workflow.basedir}/envs/qualimap.yml"
    shell:
        """
        mkdir -p $(dirname {output.TXT})

        qualimap rnaseq \
            -bam {input.BAM} \
            --sequencing-protocol strand-specific-forward \
            --sorted \
            -gtf {params.GENES_GTF} \
            --java-mem-size={resources.mem} \
            -outdir $(dirname {output.TXT}) \
            -outformat html \
        2> {log.log}
        """


# QC on STAR outputs
## qualimap on deduplicated/aligned reads
rule qualimapQC_rRNA_STAR:
    input:
        BAM="{OUTDIR}/{SAMPLE}/rRNA/STARsolo/Aligned.sortedByCoord.out.bam",
        BAI="{OUTDIR}/{SAMPLE}/rRNA/STARsolo/Aligned.sortedByCoord.out.bam.bai",
    output:
        TXT="{OUTDIR}/{SAMPLE}/qualimap/rRNA/STAR/rnaseq_qc_results.txt",
        HTML="{OUTDIR}/{SAMPLE}/qualimap/rRNA/STAR/qualimapReport.html",
    params:
        GENES_GTF=lambda wildcards: SAMPLE_SHEET["rRNA_gtf"][wildcards.SAMPLE],
    resources:
        mem="32G",
        threads=1,
    conda:
        f"{workflow.basedir}/envs/qualimap.yml"
    shell:
        """
        mkdir -p $(dirname {output.TXT})

        qualimap rnaseq \
            -bam {input.BAM} \
            --sequencing-protocol strand-specific-forward \
            --sorted \
            -gtf {params.GENES_GTF} \
            --java-mem-size={resources.mem} \
            -outdir $(dirname {output.TXT}) \
            -outformat html
        """


rule qualimap_summary2csv_rRNA_STAR:
    input:
        TXT="{OUTDIR}/{SAMPLE}/qualimap/rRNA/{TOOL}/rnaseq_qc_results.txt",
    output:
        CSV="{OUTDIR}/{SAMPLE}/qualimap/rRNA/{TOOL}/rnaseq_qc_results.csv",
    resources:
        mem="4G",
        threads=1,
    shell:
        """
        python scripts/py/qualimap_summary2csv.py {input.TXT} {output.CSV}
        """
