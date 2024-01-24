# Filter .fastqs to only keep short reads (for small RNA stuff)
rule shortPass_filter_fromfastq:
    input:
        R1_FQ = '{OUTDIR}/{sample}/tmp/final_R1.fq.gz',
        R2_FQ = '{OUTDIR}/{sample}/tmp/final_R2.fq.gz',
        R1_FQ_FILTERED = '{OUTDIR}/{sample}/tmp/final_filtered_R1.fq.gz',
        R2_FQ_FILTERED = '{OUTDIR}/{sample}/tmp/final_filtered_R2.fq.gz',
    output:        
        R1_FQ_FILTERED = '{OUTDIR}/{sample}/tmp/final_filtered_short_R1.fq.gz',
        R2_FQ_FILTERED = '{OUTDIR}/{sample}/tmp/final_filtered_short_R2.fq.gz',
    params:
        OUTDIR = config['OUTDIR'],
        MEMLIMIT = config['MEMLIMIT'],
        MAX_SHORT_READ_LENGTH = config['MAX_SHORT_READ_LENGTH']
    threads:
        config['CORES']
    priority:
        42
    run:
        recipe = RECIPE_DICT[wildcards.sample]

        # Select R2 based on alignment recipe
        if "rRNA" in recipe: # Use trimmed & rRNA-filtered .fq's
            R1 = input.R1_FQ_FILTERED
            R2 = input.R2_FQ_FILTERED
        else: # just trimmed .fq's
            R1 = input.R1_FQ
            R2 = input.R2_FQ

        # Filter paired reads where the second read is less than 40 bases 
        # We use paste to combine the two fastq files line by line on stdout
        # Then use awk to get the length of the sequence of the second read (6th column), and finally print 
        # either the first or second half of the line to the respective output file.
    shell(
        f"""
        paste <(gzip -dc {R1}) <(gzip -dc {R2}) | awk '
            {{
                l = length(\$6);
                if (l < {params.MAX_SHORT_READ_LENGTH}) {{
                    print \$1"\\n"\$2"\\n"\$3"\\n"\$4 | "gzip > {output.R1_FQ_FILTERED}"
                    print \$5"\\n"\$6"\\n"\$7"\\n"\$8 | "gzip > {output.R2_FQ_FILTERED}"
                }}
            }}
        '
        """
    )
        