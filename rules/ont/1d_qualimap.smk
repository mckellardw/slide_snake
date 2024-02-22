
rule ont_qualimap:
    input:
        BAM = "{OUTDIR}/{SAMPLE}/ont/sorted.bam", 
        BAI = "{OUTDIR}/{SAMPLE}/ont/sorted.bam.bai",
    output:
        TXT  = '{OUTDIR}/{SAMPLE}/qualimap/ont/minimap2/rnaseq_qc_results.txt',
        HTML = '{OUTDIR}/{SAMPLE}/qualimap/ont/minimap2/qualimapReport.html'
    params:
        GENES_GTF = lambda wildcards: GTF_DICT[wildcards.SAMPLE]
    threads:
        1
    resources:
        mem = "32G"
    run:
        shell(
            f"""
            mkdir -p $(dirname {output.TXT})

            {EXEC['QUALIMAP']} rnaseq \
                -bam {input.BAM} \
                -gtf {params.GENES_GTF} \
                --sequencing-protocol strand-specific-forward \
                --sorted \
                --java-mem-size={resources.mem} \
                -outdir $(dirname {output.TXT}) \
                -outformat html
            """ 
        )


rule ont_qualimap_summary2csv:
    input:
        TXT = '{OUTDIR}/{SAMPLE}/qualimap/ont/{TOOL}/rnaseq_qc_results.txt',
    output:
        CSV = '{OUTDIR}/{SAMPLE}/qualimap/ont/{TOOL}/rnaseq_qc_results.csv',
    threads:
        1
    run:
        shell(
            f"""
            python scripts/py/qualimap_summary2csv.py {input.TXT} {output.CSV}
            """
        )