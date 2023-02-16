# Remove reads that don't have a corrected spot/cell barcode with samtools, then remove duplicates w/ **umi-tools**
## High mem usage? Check here! https://umi-tools.readthedocs.io/en/latest/faq.html
#TODO: split bam by strand and by chromosome, then dedup each chr!
rule umitools_dedupBAM:
    input:
        BB = "{OUTDIR}/{sample}/bb/whitelist.txt",
        SORTEDBAM = '{OUTDIR}/{sample}/STARsolo/Aligned.sortedByCoord.out.bam'
    output:
        DEDUPBAM = '{OUTDIR}/{sample}/STARsolo/Aligned.sortedByCoord.dedup.out.bam',
        TMPBAM = temp('{OUTDIR}/{sample}/tmp/tmp.bam')
    params:
        OUTPUT_PREFIX='{OUTDIR}/{sample}/umitools_dedup/{sample}',
        # TMPBAM = '{OUTDIR}/{sample}/tmp.bam'
    threads:
        config['CORES']
        #1
    conda:
        "STARsolo"
    log:
        '{OUTDIR}/{sample}/umitools_dedup/dedup.log'
    shell:
        """
        samtools view -1 -b \
        -@ {threads} \
        --tag-file CB:{input.BB} \
        {input.SORTEDBAM} \
        > {output.TMPBAM}

        samtools index \
        -@ {threads} \
        {output.TMPBAM}

        umi_tools dedup \
        -I {output.TMPBAM} \
        --extract-umi-method=tag \
        --umi-tag=UB \
        --cell-tag=CB \
        --method=unique \
        --per-cell \
        --unmapped-reads=discard \
        --output-stats={params.OUTPUT_PREFIX} \
        --log {log} \
        -S {output.DEDUPBAM}
        """
        # rm {params.TMPBAM}
        # rm (params.TMPBAM).bai

rule umitools_indexDedupBAM:
    input:
        SORTEDBAM = '{OUTDIR}/{sample}/Aligned.sortedByCoord.dedup.out.bam'
    output:
        BAI = '{OUTDIR}/{sample}/Aligned.sortedByCoord.dedup.out.bam.bai'
    threads:
        config['CORES']
    conda:
        "STARsolo"
    shell:
        """
        samtools index -@ {threads} {input.SORTEDBAM}
        """
