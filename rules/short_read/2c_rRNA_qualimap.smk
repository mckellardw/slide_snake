# Qualimap on bwa output
## qualimap on deduplicated/aligned reads
rule ilmn_2c_qualimapQC_rRNA_bwa:
    input:
        BAM="{OUTDIR}/{SAMPLE}/short_read/rRNA/bwa/aligned_sorted.bam",
        BAI="{OUTDIR}/{SAMPLE}/short_read/rRNA/bwa/aligned_sorted.bam.bai",
    output:
        TXT="{OUTDIR}/{SAMPLE}/short_read/qualimap/rRNA/bwa/rnaseq_qc_results.txt",
        HTML="{OUTDIR}/{SAMPLE}/short_read/qualimap/rRNA/bwa/report.pdf",
    params:
        GENES_GTF=lambda wildcards: SAMPLE_SHEET["rRNA_gtf"][wildcards.SAMPLE],
    log:
        log="{OUTDIR}/{SAMPLE}/short_read/qualimap/rRNA/bwa/rnaseq_qc.log",
        err="{OUTDIR}/{SAMPLE}/short_read/qualimap/rRNA/bwa/rnaseq_qc.err",
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
            --sequencing-protocol strand-specific-forward \
            --sorted \
            -gtf {params.GENES_GTF} \
            --java-mem-size={resources.mem} \
            -outdir $(dirname {output.TXT}) \
            -outformat pdf \
        1> {log.log}
        2> {log.err}
        """


# Convert qualimap summary to readable .csv
rule ilmn_2c_qualimap_summary2csv_rRNA_STAR:
    input:
        TXT="{OUTDIR}/{SAMPLE}/short_read/qualimap/rRNA/{TOOL}/rnaseq_qc_results.txt",
    output:
        CSV="{OUTDIR}/{SAMPLE}/short_read/qualimap/rRNA/{TOOL}/rnaseq_qc_results.csv",
    resources:
        mem="4G",
    threads: 1
    shell:
        """
        python scripts/py/qualimap_summary2csv.py {input.TXT} {output.CSV}
        """
