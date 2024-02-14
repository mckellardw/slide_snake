# QC on STAR outputs

## qualimap on deduplicated/aligned reads
rule qualimapQC_STAR:
    input:
        SORTEDBAM = '{OUTDIR}/{SAMPLE}/STARsolo/{RECIPE}/Aligned.sortedByCoord.dedup.out.bam'
    output:
        qualimapDir = directory('{OUTDIR}/{SAMPLE}/qualimap/{RECIPE}'),
        qualimapReport_txt = '{OUTDIR}/{SAMPLE}/qualimap/{RECIPE}/rnaseq_qc_results.txt',
        qualimapReport_html = '{OUTDIR}/{SAMPLE}/qualimap/{RECIPE}/report.html'
    params:
        GENES_GTF = lambda wildcards: GTF_DICT[wildcards.SAMPLE]
    threads:
        1
    run:
        shell(
            f"""
            mkdir -p {output.qualimapDir}

            {EXEC['QUALIMAP']} rnaseq \
                -bam {input.SORTEDBAM} \
                -gtf {params.GENES_GTF} \
                --sequencing-protocol strand-specific-forward \
                --sorted \
                --java-mem-size=8G \
                -outdir {output.qualimapDir} \
                -outformat html
            """ 
        )
        # cd {output.qualimapDir}
        # -nt {threads} \


rule qualimap_summary2csv_STAR:
    input:
        TXT = '{OUTDIR}/{SAMPLE}/qualimap/{RECIPE}/rnaseq_qc_results.txt'
    output:
        CSV = '{OUTDIR}/{SAMPLE}/qualimap/{RECIPE}/rnaseq_qc_result.csv'
    threads:
        1
    run:
        shell(
            f"""
            python scripts/py/qualimap_summary2csv.py {input.TXT} {output.CSV}
            """
        )
