# QC on STAR outputs

## qualimap on deduplicated/aligned reads
rule qualimapQC_rRNA_STAR:
    input:
        BAM = '{OUTDIR}/{SAMPLE}/rRNA/STARsolo/Aligned.sortedByCoord.out.bam', 
        BAI = '{OUTDIR}/{SAMPLE}/rRNA/STARsolo/Aligned.sortedByCoord.out.bam.bai',
    output:
        TXT  = '{OUTDIR}/{SAMPLE}/qualimap/rRNA/STAR/rnaseq_qc_results.txt',
        HTML = '{OUTDIR}/{SAMPLE}/qualimap/rRNA/STAR/report.html'
    params:
        # GENES_GTF = lambda wildcards: GTF_DICT[wildcards.SAMPLE]
        GENES_GTF = '/gpfs/commons/groups/innovation/dwm/ref_snake/out/mus_musculus/rRNA/raw/annotations.gtf' #TODO
    threads:
        1
    run:
        shell(
            f"""
            mkdir -p $(dirname {output.TXT})

            {EXEC['QUALIMAP']} rnaseq \
                -bam {input.BAM} \
                -gtf {params.GENES_GTF} \
                --sequencing-protocol strand-specific-forward \
                --sorted \
                --java-mem-size=8G \
                -outdir $(dirname {output.TXT}) \
                -outformat html
            """ 
        )
        # cd {output.qualimapDir}
        # -nt {threads} \


# Qualimap on bwa outputs
## qualimap on deduplicated/aligned reads
rule qualimapQC_rRNA_bwa:
    input:
        BAM = '{OUTDIR}/{SAMPLE}/rRNA/bwa/aligned_sorted.bam', 
        BAI = '{OUTDIR}/{SAMPLE}/rRNA/bwa/aligned_sorted.bam.bai',
    output:
        TXT = '{OUTDIR}/{SAMPLE}/qualimap/rRNA/bwa/rnaseq_qc_results.txt',
        HTML = '{OUTDIR}/{SAMPLE}/qualimap/rRNA/bwa/report.html'
    params:
        # GENES_GTF = lambda wildcards: GTF_DICT[wildcards.SAMPLE]
        GENES_GTF = '/gpfs/commons/groups/innovation/dwm/ref_snake/out/mus_musculus/rRNA/raw/annotations.gtf' #TODO
    threads:
        1
    run:
        shell(
            f"""
            mkdir -p $(dirname {output.TXT})

            {EXEC['QUALIMAP']} rnaseq \
                -bam {input.BAM} \
                -gtf {params.GENES_GTF} \
                --sequencing-protocol strand-specific-forward \
                --sorted \
                --java-mem-size=8G \
                -outdir $(dirname {output.TXT}) \
                -outformat html
            """ 
        )
        # cd {output.qualimapDir}
        # -nt {threads} \


rule qualimap_summary2csv_rRNA_STAR:
    input:
        TXT = '{OUTDIR}/{SAMPLE}/qualimap/rRNA/{TOOL}/rnaseq_qc_results.txt',
    output:
        CSV = '{OUTDIR}/{SAMPLE}/qualimap/rRNA/{TOOL}/rnaseq_qc_results.csv',
    threads:
        1
    run:
        shell(
            f"""
            python scripts/py/qualimap_summary2csv.py {input.TXT} {output.CSV}
            """
        )