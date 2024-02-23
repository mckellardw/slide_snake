### helper functions ########################################################################
def get_barcode_length(w):
    """
    Get barcode length based on the recipe(s) passed.
    """
    bc_lengths = [RECIPE_SHEET["BC.length"][R] for R in RECIPE_DICT[wildcards.SAMPLE]]
    if len(bc_lengths) == 1:
        out = bc_lengths[0]
    else:
        #TODO: there is probably a better way to handle multi-recipe than this
        out = max(set(bc_lengths), key=bc_lengths.count) # return mode
    return out


def get_umi_length(w):
    """
    Get UMI length based on the recipe(s) passed.
    """
    umi_lengths = [RECIPE_SHEET["UMI.length"][R] for R in RECIPE_DICT[wildcards.SAMPLE]]
    if len(umi_lengths) == 1:
        out = umi_lengths[0]
    else:
        #TODO: there is probably a better way to handle multi-recipe than this
        out = max(set(umi_lengths), key=umi_lengths.count) # return mode
    return out


def get_split_ont_align_mem_gb(wildcards, threads):
    return config["RESOURCES_MM2_MEM_GB"] / threads

### rules ###################################################################################

# Borrowed portions of this code from `sockeye` - https://github.com/nanoporetech/sockeye/tree/master
rule merge_formats_ONT:
    output:
        MERGED_FQ = temp("{OUTDIR}/{SAMPLE}/ont/merged.fq.gz")
    params:
        TMPDIR = "{OUTDIR}/{SAMPLE}/ont/tmp",
        ONT_reads = lambda wildcards: ONT[wildcards.SAMPLE]
    log:
        log = "{OUTDIR}/{SAMPLE}/ont/merged.log"
    threads:
        config['CORES']
    run:
        shell(f"rm -rf {params.TMPDIR}/*") # clear out tmp dir
        if len(params.ONT_reads) == 1:
            F = params.ONT_reads[0]
            if ".fq.gz" in F or ".fastq.gz" in F:
                shell(
                    f"""
                    mkdir -p {params.TMPDIR}
                    cp {F} {output.MERGED_FQ}
                    """
                )
            elif ".sam" in F or "bam" in F:
                shell(
                    f"""
                    if [ -f {F} ]; then
                        mkdir -p {params.TMPDIR}

                        {EXEC['SAMTOOLS']} fastq {F} \
                        > {params.TMPDIR}/{F_base}.fq \
                        2> {log.log}
                    else
                        echo "File [ {F} ] does not exist." >> {log.log}
                    fi                    
                    """
                )
        else:
            for F in params.ONT_reads:
                F_base = os.path.basename(F).split('.')[0]
                if ".fq.gz" in F or ".fastq.gz" in F:
                    shell(
                        f"""
                        if [ -f {F} ]; then
                            mkdir -p {params.TMPDIR}
                            echo "Adding {F} to output fastq" >> {log.log}

                            zcat {F} > {params.TMPDIR}/{F_base}.fq
                        else
                            echo "File [ {F} ] does not exist." >> {log.log}
                        fi
                        """
                    )
                elif ".sam" in F or "bam" in F:
                    shell(
                        f"""
                        if [ -f {F} ]; then
                            mkdir -p {params.TMPDIR}

                            {EXEC['SAMTOOLS']} fastq {F} > {params.TMPDIR}/{F_base}.fq \
                            2> {log.log}
                        else
                            echo "File [ {F} ] does not exist." >> {log.log}
                        fi                    
                        """
                    )
                # end loop

        shell(
            f"""
            cat {params.TMPDIR}/*.fq > {output.MERGED_FQ.strip('.gz')}
            {EXEC['PIGZ']} -p{threads} {output.MERGED_FQ.strip('.gz')}
            """
        )
    #


# borrowed from sockeye
rule ont_call_adapter_scan:
    input:
        FQ = "{OUTDIR}/{SAMPLE}/ont/merged.fq.gz"
        # FOFN = "{OUTDIR}/{SAMPLE}/ont/chunk/fofn.txt",
    output:
        TSV = "{OUTDIR}/{SAMPLE}/ont/adapter_scan.tsv",
        FQ  = "{OUTDIR}/{SAMPLE}/ont/merged_stranded.fq.gz"
    threads: 
        config["CORES"]
    params:
        # batch_size = config["READ_STRUCTURE_BATCH_SIZE"],
        KIT = "3prime" #['3prime', '5prime', 'multiome']
        # kit=lambda w: sample_sheet.loc[w.run_id, "kit_name"]
    # conda:
    #     "../envs/stranding.yml"
    run:
        shell(
            f"""
            python scripts/py/adapter_scan_vsearch.py \
                --kit {params.KIT} \
                --output_fastq {output.FQ} \
                --output_tsv {output.TSV} \
                -t {threads} \
                {input.FQ} 
            """
        )
# --batch_size {params.batch_size} \

# Write lists of read IDs for each adapter type 
rule ont_readIDs_by_adapter_type:
    input:
        TSV = "{OUTDIR}/{SAMPLE}/ont/adapter_scan.tsv",
        FQ  = "{OUTDIR}/{SAMPLE}/ont/merged_stranded.fq.gz"
    output:
        LST = "{OUTDIR}/{SAMPLE}/ont/adapter_scan_readids/full_len.txt",
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

#TODO
rule ont_split_fastq_by_adapter_type:
    input:
        FQ  = "{OUTDIR}/{SAMPLE}/ont/merged_stranded.fq.gz",
        LST = "{OUTDIR}/{SAMPLE}/ont/adapter_scan_readids/full_len.txt",
        # DIR = "{OUTDIR}/{SAMPLE}/ont/adapter_scan_readids",
    output:
        FQ  = temp("{OUTDIR}/{SAMPLE}/ont/merged_stranded.fq"),
        SUBFQ = "{OUTDIR}/{SAMPLE}/ont/adapter_scan_readids/full_len.fq.gz"
    threads: 
        config["CORES"]
    log:
        log = "{OUTDIR}/{SAMPLE}/ont/adapter_scan_readids/subseq_full_len.log"
    run:
        # for ADAPTER in input.ADAPTER_TYPES: #TODO- broaden to other read types, bneyond full_len
        shell(
            f"""
            mkdir -p $(dirname {output.FQ})
            
            zcat {input.FQ} > {output.FQ} 

            {EXEC['SEQTK']} subseq \
                {output.FQ} \
                {input.LST} \
            > {output.SUBFQ.strip('.gz')} \
            2> {log.log}
            """
        )

        shell(
            f"""
            {EXEC['PIGZ']} -p{threads} {output.SUBFQ.strip('.gz')}
            """
        )

#TODO- rule to split fastqs into R1 + R2, then input to STARsolo
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