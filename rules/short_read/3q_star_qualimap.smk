# ALignment qc with qualimap
## link: https://qualimap.conesalab.org/


## qualimap on "raw" (not deduplicated) aligned reads
rule ilmn_3q_qualimapQC_STAR:
    input:
        BAM="{OUTDIR}/{SAMPLE}/short_read/STARsolo/{RECIPE}/Aligned.sortedByCoord.out.bam",
        BAI="{OUTDIR}/{SAMPLE}/short_read/STARsolo/{RECIPE}/Aligned.sortedByCoord.out.bam.bai",
    output:
        TXT="{OUTDIR}/{SAMPLE}/short_read/qualimap/STAR/{RECIPE}/raw/rnaseq/rnaseq_qc_results.txt",
        PDF="{OUTDIR}/{SAMPLE}/short_read/qualimap/STAR/{RECIPE}/raw/rnaseq/report.pdf",
    params:
        GENES_GTF=lambda wildcards: SAMPLE_SHEET["genes_gtf"][wildcards.SAMPLE],
    log:
        log="{OUTDIR}/{SAMPLE}/short_read/qualimap/STAR/{RECIPE}/raw/rnaseq/rnaseq.log",
        err="{OUTDIR}/{SAMPLE}/short_read/qualimap/STAR/{RECIPE}/raw/rnaseq/rnaseq.err",
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
        1> {log.log}
        2> {log.err}
        """


## qualimap on deduplicated/aligned reads
rule ilmn_3q_qualimapQC_dedup_STAR:
    input:
        BAM="{OUTDIR}/{SAMPLE}/short_read/STARsolo/{RECIPE}/Aligned.sortedByCoord.out.dedup.bam",
        BAI="{OUTDIR}/{SAMPLE}/short_read/STARsolo/{RECIPE}/Aligned.sortedByCoord.out.dedup.bam.bai",
    output:
        TXT="{OUTDIR}/{SAMPLE}/short_read/qualimap/STAR/{RECIPE}/dedup/rnaseq/rnaseq_qc_results.txt",
        PDF="{OUTDIR}/{SAMPLE}/short_read/qualimap/STAR/{RECIPE}/dedup/rnaseq/report.pdf",
    params:
        GENES_GTF=lambda wildcards: SAMPLE_SHEET["genes_gtf"][wildcards.SAMPLE],
    log:
        log="{OUTDIR}/{SAMPLE}/short_read/qualimap/STAR/{RECIPE}/dedup/rnaseq/rnaseq.log",
        err="{OUTDIR}/{SAMPLE}/short_read/qualimap/STAR/{RECIPE}/dedup/rnaseq/rnaseq.err",
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
        1> {log.log}
        2> {log.err}
        """


rule ilmn_3q_qualimap_readqc_summary2csv_STAR:
    input:
        TXT="{OUTDIR}/{SAMPLE}/short_read/qualimap/STAR/{RECIPE}/{DEDUP}/rnaseq/rnaseq_qc_results.txt",
    output:
        CSV="{OUTDIR}/{SAMPLE}/short_read/qualimap/STAR/{RECIPE}/{DEDUP}/rnaseq/rnaseq_qc_results.csv",
    threads: 1
    shell:
        """
        python scripts/py/qualimap_summary2csv.py \
            --input {input.TXT} \
            --output {output.CSV}
        """


rule ilmn_3q_qualimap_bamqc_STAR_raw:
    input:
        BAM="{OUTDIR}/{SAMPLE}/short_read/STARsolo/{RECIPE}/Aligned.sortedByCoord.out.bam",
        BAI="{OUTDIR}/{SAMPLE}/short_read/STARsolo/{RECIPE}/Aligned.sortedByCoord.out.bam.bai",
    output:
        REPORT="{OUTDIR}/{SAMPLE}/short_read/qualimap/STAR/{RECIPE}/raw/bamqc/report.pdf",
        CSV="{OUTDIR}/{SAMPLE}/short_read/qualimap/STAR/{RECIPE}/raw/bamqc/genome_results.txt",
    log:
        log="{OUTDIR}/{SAMPLE}/short_read/qualimap/STAR/{RECIPE}/raw/bamqc/bamqc.log",
        err="{OUTDIR}/{SAMPLE}/short_read/qualimap/STAR/{RECIPE}/raw/bamqc/bamqc.err",
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


rule ilmn_3q_qualimap_bamqc_STAR_dedup:
    input:
        BAM="{OUTDIR}/{SAMPLE}/short_read/STARsolo/{RECIPE}/Aligned.sortedByCoord.out.dedup.bam",
        BAI="{OUTDIR}/{SAMPLE}/short_read/STARsolo/{RECIPE}/Aligned.sortedByCoord.out.dedup.bam.bai",
    output:
        REPORT="{OUTDIR}/{SAMPLE}/short_read/qualimap/STAR/{RECIPE}/dedup/bamqc/report.pdf",
        CSV="{OUTDIR}/{SAMPLE}/short_read/qualimap/STAR/{RECIPE}/dedup/bamqc/genome_results.txt",
    log:
        log="{OUTDIR}/{SAMPLE}/short_read/qualimap/STAR/{RECIPE}/dedup/bamqc/bamqc.log",
        err="{OUTDIR}/{SAMPLE}/short_read/qualimap/STAR/{RECIPE}/dedup/bamqc/bamqc.err",
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


rule ilmn_3q_qualimap_bamqc_summary2csv_STAR:
    input:
        TXT="{OUTDIR}/{SAMPLE}/short_read/qualimap/STAR/{RECIPE}/{DEDUP}/bamqc/genome_results.txt",
    output:
        CSV="{OUTDIR}/{SAMPLE}/short_read/qualimap/STAR/{RECIPE}/{DEDUP}/bamqc/genome_results.csv",
    threads: 1
    shell:
        """
        python scripts/py/qualimap_summary2csv.py \
            --input {input.TXT} \
            --output {output.CSV}
        """
