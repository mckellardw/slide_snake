# fastqc before trimming
rule ont_readQC_preCutadapt:
    input:
        FQ="{OUTDIR}/{SAMPLE}/tmp/ont/adapter_scan_readids/merged_adapter_{READ}.fq.gz",
    output:
        TSV="{OUTDIR}/{SAMPLE}/ont/readqc/1_preCutadapt/{READ}_qc.tsv", #TODO compress this?
    log:
        log="{OUTDIR}/{SAMPLE}/ont/readqc/1_preCutadapt/{READ}_qc.log"
    threads: config["CORES"]
    shell:
        """
        python scripts/py/fastq_read_qc.py \
            {input.FQ} \
            {output.TSV} \
            --threads {threads} \
        2>&1 | tee {log.log}
        """


# fastqc after cutadapt trimming
rule ont_readQC_postCutadapt:
    input:
        FQ="{OUTDIR}/{SAMPLE}/tmp/ont/cut_{READ}.fq.gz",
    output:
        TSV="{OUTDIR}/{SAMPLE}/ont/readqc/2_postCutadapt/{READ}_qc.tsv",
    log:
        log="{OUTDIR}/{SAMPLE}/ont/readqc/2_postCutadapt/{READ}_qc.log"
    threads: config["CORES"]
    shell:
        """
        python scripts/py/fastq_read_qc.py \
            {input.FQ} \
            {output.TSV} \
            --threads {threads} \
        2>&1 | tee {log.log}
        """

rule ont_readQC_summaryplot:
    input:
        TSV="{OUTDIR}/{SAMPLE}/ont/readqc/{CUT}/{READ}_qc.tsv",
    output:
        IMG="{OUTDIR}/{SAMPLE}/ont/readqc/{CUT}/{READ}_qc.png",
    threads: 1    
    conda:
        f"{workflow.basedir}/envs/ggplot2.yml"
    shell:
        """
        Rscript scripts/R/readqc_summary.R -f {input.TSV} -o {output.IMG}
        """