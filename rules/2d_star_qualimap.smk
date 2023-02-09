
#############################################
## QC on STAR outputs
#############################################

## qualimap on aligned reads
#TODO- switch to dedup'ed .bam, once I re-write umitools deduplication
rule qualimapQC:
    input:
        SORTEDBAM = '{OUTDIR}/{sample}/STARsolo/Aligned.sortedByCoord.out.bam'
    output:
        qualimapDir = directory('{OUTDIR}/{sample}/qualimap'),
        fastqcReport = '{OUTDIR}/{sample}/qualimap/qualimapReport.html'
    params:
        GENES_GTF = lambda wildcards: GTF_DICT[wildcards.sample]
    threads:
        1
    conda:
        "STARsolo"
    shell:
        """
        mkdir -p {output.qualimapDir}
        cd {output.qualimapDir}

        qualimap rnaseq \
        -bam {input.SORTEDBAM} \
        -gtf {params.GENES_GTF} \
        --sequencing-protocol strand-specific-forward \
        --sorted \
        --java-mem-size=8G \
        -outdir {output.qualimapDir} \
        -outformat html
        """
        # -nt {threads} \
