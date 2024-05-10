# Borrowed portions of this code from `sockeye` - https://github.com/nanoporetech/sockeye/tree/master
# TODO- rewrite as a python script...
rule merge_formats_ONT:
    output:
        MERGED_FQ=temp("{OUTDIR}/{SAMPLE}/tmp/ont/merged.fq.gz"),
    params:
        TMPDIR="{OUTDIR}/{SAMPLE}/tmp/ont",
        ONT_reads=lambda wildcards: ONT[wildcards.SAMPLE],
        CHUNK_SIZE=50,
    log:
        log="{OUTDIR}/{SAMPLE}/ont/misc_logs/merged.log",
    threads: config["CORES"]
    run:
        shell(f"mkdir -p {params.TMPDIR}")
        # shell(f"rm -rf {params.TMPDIR}/*") # clear out tmp dir
        shell(f"echo 'Read files:' > {log.log}")
        for f in params.ONT_reads:
            shell(f"echo '   {f}' >> {log.log}")

        if len(params.ONT_reads) == 1 and "*" not in params.ONT_reads[0]:
            F = params.ONT_reads[0]
            if ".fq.gz" in F or ".fastq.gz" in F:
                shell(
                    f"""
                    cp {F} {output.MERGED_FQ} 2>> {log.log}
                    """
                )
            elif ".sam" in F or ".bam" in F:
                shell(
                    f"""
                    if [ -f {F} ]; then
                        {EXEC['SAMTOOLS']} fastq {F} \
                        > {output.MERGED_FQ.strip('.gz')} \
                        2>> {log.log}
                    else
                        echo "File [ {F} ] does not exist." >> {log.log}
                    fi                    
                    
                    {EXEC['PIGZ']} -p{threads} {output.MERGED_FQ.strip('.gz')} 2>> {log.log}
                    """
                )
        elif len(params.ONT_reads) == 1 and "*" in params.ONT_reads[0]:
            import glob

            F_list = glob.glob(params.ONT_reads[0])

            shell(f"echo 'Regex-ed file list:' >> {log.log}")
            for f in F_list:
                shell(f"echo '   {f}' >> {log.log}")

                # TODO- fix this so bam files can be passed individually
                # if len(F_list) > params.CHUNK_SIZE:
                #     shell(f"echo 'Chunking regex-ed file list...' >> {log.log}")
                #     F_list_chunked = [
                #         " ".join(F_list[i:i + params.CHUNK_SIZE])
                #             for i in range(0, len(F_list),  params.CHUNK_SIZE)
                #     ]
                # else:
                #     F_list_chunked = " ".join(F_list)

            F_base = output.MERGED_FQ.strip(".gz")

            for i, F in enumerate(F_list):  # F_list_chunked
                if ".fq.gz" in F or ".fastq.gz" in F:
                    shell(
                        f"""
                        echo "Adding {F} to output fastq" >> {log.log}

                        zcat {F} >> {F_base}
                        """
                    )
                elif ".sam" in F or ".bam" in F:
                    shell(
                        f"""
                        {EXEC['SAMTOOLS']} fastq {F} \
                        >> {F_base} \
                        2>> {log.log}
                        """
                    )
                    # end for loop

            shell(
                f"""
                {EXEC['PIGZ']} -p{threads} {F_base}
                """
            )

        else:
            F_base = output.MERGED_FQ.strip(".gz")
            for F in params.ONT_reads:
                if ".fq.gz" in F or ".fastq.gz" in F:
                    shell(
                        f"""
                        if [ -f {F} ]; then
                            echo "Adding {F} to output fastq" >> {log.log}

                            zcat {F} >> {F_base}
                        else
                            echo "File [ {F} ] does not exist." >> {log.log}
                        fi
                        """
                    )
                elif ".sam" in F or ".bam" in F:
                    shell(
                        f"""
                        if [ -f {F} ]; then
                            {EXEC['SAMTOOLS']} fastq {F} \
                            >> {F_base} \
                            2>> {log.log}
                        else
                            echo "File [ {F} ] does not exist." >> {log.log}
                        fi                    
                        """
                    )
                    # end loop

            shell(
                f"""
                {EXEC['PIGZ']} -p{threads} {F_base}  2>> {log.log}
                """
            )


# borrowed from sockeye
rule ont_call_adapter_scan:
    input:
        FQ="{OUTDIR}/{SAMPLE}/tmp/ont/merged.fq.gz",
    output:
        TSV="{OUTDIR}/{SAMPLE}/ont/adapter_scan.tsv",
        FQ="{OUTDIR}/{SAMPLE}/tmp/ont/merged_stranded.fq.gz",
    threads: config["CORES"]
    params:
        # batch_size = config["READ_STRUCTURE_BATCH_SIZE"],
        KIT="3prime",  #['3prime', '5prime', 'multiome']
    # conda:
    #     f"{workflow.basedir}/envs/slsn_ont_prep.yml"
    log:
        log="{OUTDIR}/{SAMPLE}/ont/misc_logs/adapter_scan.log",
    shell:
        """
        python scripts/py/adapter_scan_vsearch.py \
            --kit {params.KIT} \
            --output_fastq {output.FQ} \
            --output_tsv {output.TSV} \
            -t {threads} \
            {input.FQ} \
        2>&1 | tee {log.log}
        """


# Write lists of read IDs for each adapter type
rule ont_readIDs_by_adapter_type:
    input:
        TSV="{OUTDIR}/{SAMPLE}/ont/adapter_scan.tsv",
        FQ="{OUTDIR}/{SAMPLE}/tmp/ont/merged_stranded.fq.gz",
    output:
        FULL_LEN="{OUTDIR}/{SAMPLE}/ont/adapter_scan_readids/full_len.txt",
        SINGLE_ADAPTER1="{OUTDIR}/{SAMPLE}/ont/adapter_scan_readids/single_adapter1.txt",
    threads: config["CORES"]
    run:
        shell(
            f"""
            python scripts/py/write_adapterscan_read_id_lists.py \
                --tsv_file_path {input.TSV} \
                --output_directory $(dirname {output.FULL_LEN})
            """
        )


# Write lists of read IDs for each adapter type
rule ont_adapterScan_QC:
    input:
        FULL_LEN="{OUTDIR}/{SAMPLE}/ont/adapter_scan_readids/full_len.txt",
        SINGLE_ADAPTER1="{OUTDIR}/{SAMPLE}/ont/adapter_scan_readids/single_adapter1.txt",
    output:
        LOG="{OUTDIR}/{SAMPLE}/ont/misc_logs/adapter_scan_results.txt",
    threads: 1
    run:
        shell(
            f"""
            dir_path=$(dirname {input.FULL_LEN})

            for file in "$dir_path"/*.txt; do
                echo "$(basename $file)"\t"$(wc -l <"$file")" >> {output.LOG}
            done
            """
        )


# merge lists of useable reads
## FULL_LEN = R1 sequence & TSO sequence
## SINGLE_ADAPTER1 = just R1 sequence - incompletely sequenced
rule ont_merge_scan_lists:
    input:
        LOG="{OUTDIR}/{SAMPLE}/ont/misc_logs/adapter_scan_results.txt",
        FULL_LEN="{OUTDIR}/{SAMPLE}/ont/adapter_scan_readids/full_len.txt",
        SINGLE_ADAPTER1="{OUTDIR}/{SAMPLE}/ont/adapter_scan_readids/single_adapter1.txt",
    output:
        LST="{OUTDIR}/{SAMPLE}/ont/adapter_scan_readids/keep.txt",
    threads: 1
    shell:
        """
        cat {input.FULL_LEN} {input.SINGLE_ADAPTER1} > {output.LST}
        """


# TODO- add more functionality for other read/adapter types to salvage imperfect reads
rule ont_subset_fastq_by_adapter_type:
    input:
        FQ="{OUTDIR}/{SAMPLE}/tmp/ont/merged_stranded.fq.gz",
        LST="{OUTDIR}/{SAMPLE}/ont/adapter_scan_readids/keep.txt",
        # FULL_LEN = "{OUTDIR}/{SAMPLE}/ont/adapter_scan_readids/full_len.txt",
        # SINGLE_ADAPTER1 = "{OUTDIR}/{SAMPLE}/ont/adapter_scan_readids/single_adapter1.txt",
    output:
        FQ=temp("{OUTDIR}/{SAMPLE}/tmp/ont/merged_stranded.fq"),
        # FULL_LEN = "{OUTDIR}/{SAMPLE}/tmp/ont/adapter_scan_readids/full_len.fq.gz",
        # SINGLE_ADAPTER1 = "{OUTDIR}/{SAMPLE}/tmp/ont/adapter_scan_readids/single_adapter1.fq.gz",
        FQ_ADAPTER="{OUTDIR}/{SAMPLE}/tmp/ont/adapter_scan_readids/merged_adapter.fq.gz",
    threads: config["CORES"]
    log:
        log="{OUTDIR}/{SAMPLE}/ont/misc_logs/subseq_full_len.log",
    run:
        # for ADAPTER in input.ADAPTER_TYPES:
        shell(
            f"""
            mkdir -p $(dirname {output.FQ})            
            zcat {input.FQ} > {output.FQ} 

            {EXEC['SEQTK']} subseq \
                {output.FQ} \
                {input.LST} \
            > {output.FQ_ADAPTER.strip('.gz')} \
            2> {log.log}
                        
            {EXEC['PIGZ']} -p{threads} {output.FQ_ADAPTER.strip('.gz')}
            """
        )


# TODO- parallelize by chunking input fastq
rule ont_split_fastq_to_R1_R2:
    input:
        FQ="{OUTDIR}/{SAMPLE}/tmp/ont/adapter_scan_readids/merged_adapter.fq.gz",
    output:
        R1_FQ="{OUTDIR}/{SAMPLE}/tmp/ont/adapter_scan_readids/merged_adapter_R1.fq.gz",
        R2_FQ="{OUTDIR}/{SAMPLE}/tmp/ont/adapter_scan_readids/merged_adapter_R2.fq.gz",
    params:
        ADAPTER="T" * 10,
    threads: config["CORES"]
    # 1
    log:
        log="{OUTDIR}/{SAMPLE}/ont/misc_logs/read_split.log",
    run:
        # for ADAPTER in input.ADAPTER_TYPES: #TODO- broaden to other read types, bneyond full_len
        shell(
            f"""
            python scripts/py/split_reads_parallelized.py \
                --fq_in {input.FQ} \
                --split_seq {params.ADAPTER} \
                --threads {threads} \
             2> {log.log}
            """
        )
