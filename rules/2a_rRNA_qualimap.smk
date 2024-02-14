# QC on STAR outputs

## qualimap on deduplicated/aligned reads
rule qualimapQC_rRNA_STAR:
    input:
        BAM = '{OUTDIR}/{SAMPLE}/rRNA/STAR/aligned_sorted.bam', 
        BAI = '{OUTDIR}/{SAMPLE}/rRNA/STAR/aligned_sorted.bam.bai',
    output:
        qualimapDir = directory('{OUTDIR}/{SAMPLE}/qualimap/{RECIPE}/rRNA/STAR/'),
        TXT = '{OUTDIR}/{SAMPLE}/qualimap/{RECIPE}/rRNA/STAR/rnaseq_qc_results.txt',
        qualimapReport_html = '{OUTDIR}/{SAMPLE}/qualimap/{RECIPE}/rRNA/STAR/qualimapReport.html'
    params:
        # GENES_GTF = lambda wildcards: GTF_DICT[wildcards.SAMPLE]
        GENES_GTF = ''
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


rule qualimap_summary2csv_rRNA_STAR:
    input:
        TXT = '{OUTDIR}/{SAMPLE}/qualimap/{RECIPE}/rnaseq_qc_results.txt'
    output:
        CSV = '{OUTDIR}/{SAMPLE}/qualimap/{RECIPE}/rnaseq_qc_result.csv'
    threads:
        1
    run:
        with open(input.TXT, 'r') as file:
            lines = file.readlines()

        out_dict = {}
        for line in lines:
            if '=' in line:
                key, value = line.split('=')
                if '(' in value:
                    value, tmp = value.split('(')
                    value = float(value.strip().replace(',',''))
                elif '%' in value:
                    value = float(value.rstrip('%').strip().replace(',',''))
                elif ',' in value:
                    value = float(value.strip().replace(',',''))
            #TODO - "Junction analysis" section
            # if ':' in line:
            #     key, value = line.split(':')
                out_dict[key.strip()]=value

        # df = pd.json_normalize(sections)
        # print(out_dict)
        out_df = pd.DataFrame.from_dict(
            out_dict, 
            orient='index'
        ) 
        out_df.T.to_csv(
            output.CSV, 
            index=False
        )

# Qualimap on bwa outputs
## qualimap on deduplicated/aligned reads
rule qualimapQC_rRNA_bwa:
    input:
        BAM = '{OUTDIR}/{SAMPLE}/rRNA/bwa/aligned_sorted.bam', 
        BAI = '{OUTDIR}/{SAMPLE}/rRNA/bwa/aligned_sorted.bam.bai',
    output:
        qualimapDir = directory('{OUTDIR}/{SAMPLE}/qualimap/{RECIPE}/rRNA/bwa/'),
        TXT = '{OUTDIR}/{SAMPLE}/qualimap/{RECIPE}/rRNA/bwa/rnaseq_qc_results.txt',
        qualimapReport_html = '{OUTDIR}/{SAMPLE}/qualimap/{RECIPE}/rRNA/bwa/qualimapReport.html'
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


rule qualimap_summary2csv_rRNA_bwa:
    input:
        TXT = '{OUTDIR}/{SAMPLE}/qualimap/{RECIPE}/rnaseq_qc_results.txt'
    output:
        CSV = '{OUTDIR}/{SAMPLE}/qualimap/{RECIPE}/rnaseq_qc_result.csv'
    threads:
        1
    run:
        with open(input.TXT, 'r') as file:
            lines = file.readlines()

        out_dict = {}
        for line in lines:
            if '=' in line:
                key, value = line.split('=')
                if '(' in value:
                    value, tmp = value.split('(')
                    value = float(value.strip().replace(',',''))
                elif '%' in value:
                    value = float(value.rstrip('%').strip().replace(',',''))
                elif ',' in value:
                    value = float(value.strip().replace(',',''))
            #TODO - "Junction analysis" section
            # if ':' in line:
            #     key, value = line.split(':')
                out_dict[key.strip()]=value

        # df = pd.json_normalize(sections)
        # print(out_dict)
        out_df = pd.DataFrame.from_dict(
            out_dict, 
            orient='index'
        ) 
        out_df.T.to_csv(
            output.CSV, 
            index=False
        )
