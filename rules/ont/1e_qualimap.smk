
rule ont_qualimap_minimap2:
    input:
        BAM = "{OUTDIR}/{SAMPLE}/ont/minimap2/{RECIPE}/sorted_bc.bam",
        BAI = "{OUTDIR}/{SAMPLE}/ont/minimap2/{RECIPE}/sorted_bc.bam.bai",
    output:
        TXT  = "{OUTDIR}/{SAMPLE}/qualimap/ont/minimap2/{RECIPE}/rnaseq_qc_results.txt",
        HTML = "{OUTDIR}/{SAMPLE}/qualimap/ont/minimap2/{RECIPE}/qualimapReport.html",
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

rule ont_qualimap_STAR:
    input:
        BAM = "{OUTDIR}/{SAMPLE}/STARsolo/ont/{RECIPE}/Aligned.sortedByCoord.out.bam", 
        BAI = "{OUTDIR}/{SAMPLE}/STARsolo/ont/{RECIPE}/Aligned.sortedByCoord.out.bam.bai"
    output:
        TXT  = "{OUTDIR}/{SAMPLE}/qualimap/ont/STARsolo/{RECIPE}/rnaseq_qc_results.txt",
        HTML = "{OUTDIR}/{SAMPLE}/qualimap/ont/STARsolo/{RECIPE}/qualimapReport.html"
    params:
        # GENES_GTF = lambda wildcards: GTF_DICT[wildcards.SAMPLE]
        GENES_GTF = "/gpfs/commons/groups/innovation/dwm/ref_snake/out/mus_musculus/genome/raw/annotations.gtf"
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
        TXT = "{OUTDIR}/{SAMPLE}/qualimap/ont/{TOOL}/rnaseq_qc_results.txt",
    output:
        CSV = "{OUTDIR}/{SAMPLE}/qualimap/ont/{TOOL}/rnaseq_qc_results.csv",
    threads:
        1
    run:
        shell(
            f"""
            python scripts/py/qualimap_summary2csv.py {input.TXT} {output.CSV}
            """
        )