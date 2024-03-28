
# Borrowed portions of this code from `sockeye` - https://github.com/nanoporetech/sockeye/tree/master
#TODO- rewrite as a python script...
rule merge_formats_ONT:
    output:
        MERGED_FQ = temp("{OUTDIR}/{SAMPLE}/tmp/ont/merged.fq.gz")
    params:
        TMPDIR = "{OUTDIR}/{SAMPLE}/tmp/ont",
        ONT_reads = lambda wildcards: ONT[wildcards.SAMPLE],
        CHUNK_SIZE = 50
    log:
        log = "{OUTDIR}/{SAMPLE}/ont/merged.log"
    threads:
        config['CORES']
    run:
        shell(f"mkdir -p {params.TMPDIR}")
        # shell(f"rm -rf {params.TMPDIR}/*") # clear out tmp dir

        if len(params.ONT_reads) == 1 and "*" not in params.ONT_reads[0]:
            F = params.ONT_reads[0]
            if ".fq.gz" in F or ".fastq.gz" in F:
                shell(
                    f"""
                    cp {F} {output.MERGED_FQ}
                    """
                )
            elif ".sam" in F or ".bam" in F:
                shell(
                    f"""
                    if [ -f {F} ]; then
                        {EXEC['SAMTOOLS']} fastq {F} \
                        > {output.MERGED_FQ.strip('.gz')} \
                        2> {log.log}
                    else
                        echo "File [ {F} ] does not exist." >> {log.log}
                    fi                    
                    
                    {EXEC['PIGZ']} -p{threads} {output.MERGED_FQ.strip('.gz')}
                    """
                )
        elif len(params.ONT_reads) == 1 and "*" in params.ONT_reads[0]:
            import glob
            
            F_list = glob.glob(params.ONT_reads[0])
            if len(F_list) > params.CHUNK_SIZE:
                F_list_chunked = [
                    " ".join(F_list[i:i + params.CHUNK_SIZE]) 
                        for i in range(0, len(F_list),  params.CHUNK_SIZE)
                ]
            else:
                F_list_chunked = " ".join(F_list)
            
            # F_base = os.path.basename(F).split('.')[0]
            # F_base.replace("*","42")
            # F_base = "tmp_star"
            F_base = output.MERGED_FQ.strip('.gz')

            for i, F in enumerate(F_list_chunked):
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
                        2> {log.log}
                        """
                )
                # end for loop

            shell(
                f"""
                {EXEC['PIGZ']} -p{threads} {F_base}
                """
            )

        else:
            F_base = output.MERGED_FQ.strip('.gz')
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
                {EXEC['PIGZ']} -p{threads} {F_base}
                """
            )


# borrowed from sockeye
rule ont_call_adapter_scan:
    input:
        FQ = "{OUTDIR}/{SAMPLE}/tmp/ont/merged.fq.gz"
        # FOFN = "{OUTDIR}/{SAMPLE}/ont/chunk/fofn.txt",
    output:
        TSV = "{OUTDIR}/{SAMPLE}/ont/adapter_scan.tsv",
        FQ  = "{OUTDIR}/{SAMPLE}/tmp/ont/merged_stranded.fq.gz"
    threads: 
        config["CORES"]
    params:
        # batch_size = config["READ_STRUCTURE_BATCH_SIZE"],
        KIT = "3prime" #['3prime', '5prime', 'multiome']
        # kit=lambda w: sample_sheet.loc[w.run_id, "kit_name"]
    # conda:
    #     f"{workflow.basedir}/envs/slsn_ont_prep.yml"
    # log:
    #     log="{OUTDIR}/{SAMPLE}/ont/adapter_scan.log"
    shell:
        """
        python scripts/py/adapter_scan_vsearch.py \
            --kit {params.KIT} \
            --output_fastq {output.FQ} \
            --output_tsv {output.TSV} \
            -t {threads} \
            {input.FQ} 
        """
# --batch_size {params.batch_size} \


# Write lists of read IDs for each adapter type 
rule ont_readIDs_by_adapter_type:
    input:
        TSV = "{OUTDIR}/{SAMPLE}/ont/adapter_scan.tsv",
        FQ  = "{OUTDIR}/{SAMPLE}/tmp/ont/merged_stranded.fq.gz"
    output:
        FULL_LEN = "{OUTDIR}/{SAMPLE}/ont/adapter_scan_readids/full_len.txt",
        SINGLE_ADAPTER1 = "{OUTDIR}/{SAMPLE}/ont/adapter_scan_readids/single_adapter1.txt",
        # DIR = directory("{OUTDIR}/{SAMPLE}/ont/adapter_scan_readids"),
    threads: 
        config["CORES"]
    run:
        shell(
            f"""
            python scripts/py/write_adapterscan_read_id_lists.py \
                --tsv_file_path {input.TSV} \
                --output_directory $(dirname {output.LST})
            """
        )

# merge lists of useable reads
## FULL_LEN = R1 sequence & TSO sequence
## SINGLE_ADAPTER1 = just R1 sequence - incompletely sequenced
rule ont_merge_scan_lists:
    input:
        FULL_LEN = "{OUTDIR}/{SAMPLE}/ont/adapter_scan_readids/full_len.txt",
        SINGLE_ADAPTER1 = "{OUTDIR}/{SAMPLE}/ont/adapter_scan_readids/single_adapter1.txt",
        # DIR = "{OUTDIR}/{SAMPLE}/ont/adapter_scan_readids",
    output:
        LST = "{OUTDIR}/{SAMPLE}/ont/adapter_scan_readids/keep.txt",
    threads: 
        1
    shell:
        """
        cat {input.FULL_LEN} {input.SINGLE_ADAPTER1} > {output.LST}
        """

#TODO- add more functionality for other read/adapter types to salvage imperfect reads
rule ont_subset_fastq_by_adapter_type:
    input:
        FQ  = "{OUTDIR}/{SAMPLE}/tmp/ont/merged_stranded.fq.gz",
        LST = "{OUTDIR}/{SAMPLE}/ont/adapter_scan_readids/keep.txt",
        # FULL_LEN = "{OUTDIR}/{SAMPLE}/ont/adapter_scan_readids/full_len.txt",
        # SINGLE_ADAPTER1 = "{OUTDIR}/{SAMPLE}/ont/adapter_scan_readids/single_adapter1.txt",
        # DIR = "{OUTDIR}/{SAMPLE}/ont/adapter_scan_readids",
    output:
        FQ  = temp("{OUTDIR}/{SAMPLE}/tmp/ont/merged_stranded.fq"),
        # FULL_LEN = "{OUTDIR}/{SAMPLE}/tmp/ont/adapter_scan_readids/full_len.fq.gz",
        # SINGLE_ADAPTER1 = "{OUTDIR}/{SAMPLE}/tmp/ont/adapter_scan_readids/single_adapter1.fq.gz",
        FQ_ADAPTER = "{OUTDIR}/{SAMPLE}/tmp/ont/adapter_scan_readids/merged_adapter.fq.gz",
    threads: 
        config["CORES"]
    log:
        log = "{OUTDIR}/{SAMPLE}/ont/adapter_scan_readids/subseq_full_len.log"
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
            # {EXEC['SEQTK']} subseq \
            #     {output.FQ} \
            #     {input.SINGLE_ADAPTER1} \
            # >> {output.FQ_ADAPTER.strip('.gz')} \
            # 2>> {log.log}

#TODO- parallelize by chunking input fastq
rule ont_split_fastq_to_R1_R2:
    input:
        FQ = "{OUTDIR}/{SAMPLE}/tmp/ont/adapter_scan_readids/merged_adapter.fq.gz",
    output:
        R1_FQ = "{OUTDIR}/{SAMPLE}/tmp/ont/adapter_scan_readids/merged_adapter_R1.fq.gz",
        R2_FQ = "{OUTDIR}/{SAMPLE}/tmp/ont/adapter_scan_readids/merged_adapter_R2.fq.gz"
    params:
        ADAPTER = "T"*10
    threads: 
        config["CORES"]
        # 1
    log:
        log = "{OUTDIR}/{SAMPLE}/ont/adapter_scan_readids/read_split.log"
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
            # {EXEC['PIGZ']} \
            #     -p{threads} \
            #     {output.R1_FQ.strip('.gz')} \
            #     {output.R2_FQ.strip('.gz')}

#TODO other rules to salvage non-full_len reads

# rule call_paftools:
#     input:
#         gtf=str(REF_GENES_GTF),
#     output:
#         bed=REF_GENES_BED,
#     # conda:
#     #     "../envs/minimap2.yml"
#     run:
#         shell(
#             f"""
#             {EXEC['K8']} scripts/js/paftools.js gff2bed \
#                 -j {input.gtf} \
#                 > {output.bed}
#             """
#         )





### TODO ################################################################################################
#TODO: Chunk fastqs for speed up
# checkpoint chunk_fastqs_ONT:
#     input:
#         MERGED_FQ = "{OUTDIR}/{SAMPLE}/ont/merged.fq.gz"#temp()
#     output:
#         FOFN = "{OUTDIR}/{SAMPLE}/ont/chunk/fofn.txt",
#         DIR = directory("{OUTDIR}/{SAMPLE}/ont/chunk")
#     params:
#         N_CHUNKS = 16
#     threads: 
#         config["CORES"]
#     # conda:
#     #     "../envs/stranding.yml"
#     run:
#         shell(
#             f"""
#             mkdir -p {output.DIR}

#             find -L {input.dir} -type f -name '*' > {output.FOFN}
#             """
#         )
        
#         shell(
#             f"""
#             python scripts/py/chunk_fastqs.py \
#                 --threads {threads} \
#                 --output_dir {output.DIR} \
#                 {input}
#             """
#         )

#TODO - add with chunking
# def gather_tsv_files_from_run(wildcards):
#     # throw and Exception if checkpoint is pending
#     checkpoint_dir = checkpoints.call_cat_fastq.get(**wildcards).output[0]
#     return expand(
#         READ_CONFIG_CHUNKED,
#         run_id=wildcards.run_id,
#         batch_id=glob_wildcards(
#             os.path.join(checkpoint_dir, "{batch_id}.fastq")
#         ).batch_id,
#     )
# rule combine_adapter_tables:
#     input:
#         tsv_files=gather_tsv_files_from_run,
#     output:
#         READ_CONFIG,
#     run:
#         import os
#         import pandas as pd

#         dfs = [pd.read_csv(fn, sep="\t") for fn in input.tsv_files]
#         df = pd.concat(dfs, axis=0)
#         df.to_csv(output[0], sep="\t", index=False)
#         [os.remove(fn) for fn in input]


# def gather_fastq_files_from_run(wildcards):
#     checkpoint_dir = checkpoints.call_cat_fastq.get(**wildcards).output[0]
#     return expand(
#         STRANDED_FQ_CHUNKED,
#         run_id=wildcards.run_id,
#         batch_id=glob_wildcards(
#             os.path.join(checkpoint_dir, "{batch_id}.fastq")
#         ).batch_id,
#     )


# rule combine_stranded_fastqs:
#     input:
#         fastq_files = gather_fastq_files_from_run,
#     output:
#         STRANDED_FQ,
#     run:
#         import shutil

#         with open(output[0], "wb") as f_out:
#             for chunked_fastq in input.fastq_files:
#                 with open(chunked_fastq, "rb") as f_:
#                     shutil.copyfileobj(f_, f_out)
#         [os.remove(fn) for fn in input]