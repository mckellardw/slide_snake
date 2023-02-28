
#############################################
## QC on STAR outputs
#############################################

## qualimap on deduplicated/aligned reads
rule qualimapQC:
    input:
        SORTEDBAM = '{OUTDIR}/{sample}/STARsolo/Aligned.sortedByCoord.dedup.out.bam'
    output:
        qualimapDir = directory('{OUTDIR}/{sample}/qualimap'),
        fastqcReport = '{OUTDIR}/{sample}/qualimap/qualimapReport.html'
    params:
        GENES_GTF = lambda wildcards: GTF_DICT[wildcards.sample]
    threads:
        1
    run:
        shell(
            f"""
            mkdir -p {output.qualimapDir}
            cd {output.qualimapDir}

            {QUALIMAP_EXEC} rnaseq \
            -bam {input.SORTEDBAM} \
            -gtf {params.GENES_GTF} \
            --sequencing-protocol strand-specific-forward \
            --sorted \
            --java-mem-size=8G \
            -outdir {output.qualimapDir} \
            -outformat html
            """ 
        )
        # -nt {threads} \
