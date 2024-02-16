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
rule call_adapter_scan:
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


def get_split_ont_align_mem_gb(wildcards, threads):
    return config["RESOURCES_MM2_MEM_GB"] / threads

rule align_ONT_minimap2:
    input:
        FQ = "{OUTDIR}/{SAMPLE}/ont/merged_stranded.fq.gz",
    output:
        SAM_TMP=temp("{OUTDIR}/{SAMPLE}/ont/tmp.sam"),
    params:
        ref = config["REF_GENOME_FASTA"],
        chrom_sizes = config["REF_CHROM_SIZES"],
        bed = config["REF_GENES_BED"],
        flags = config["RESOURCES_MM2_FLAGS"],
    log:
        log = "{OUTDIR}/{SAMPLE}/ont/minimap2.log"
    threads: 
        config["CORES"]
    # resources:
    #     mem_gb=get_split_ont_align_mem_gb,
    # conda:
    #     "../envs/minimap2.yml"
    run:
        shell(
            f"""
            {EXEC['MINIMAP2']} \
                -ax splice \
                -uf --MD \
                -t {threads} \
                --junc-bed {params.bed} \
                --secondary=no \
                {params.flags} \
                {params.ref} \
                {input.FQ} \
            2> {log.log} \
            > {output.SAM_TMP}
            """
        )
#


rule sort_index_ONT_output:
    input:
        SAM="{OUTDIR}/{SAMPLE}/ont/tmp.sam"
    output:
        # BAM_UNSORT_TMP=temp("{OUTDIR}/{SAMPLE}/ont/tmp_unsort.sam"),
        BAM="{OUTDIR}/{SAMPLE}/ont/sorted.bam",
        BAI="{OUTDIR}/{SAMPLE}/ont/sorted.bam.bai",
    params:
        ref = config["REF_GENOME_FASTA"]
    # threads: 
    #     config["CORES"]
    run:
        shell(
            f"""
            {EXEC['SAMTOOLS']} sort \
                --reference {params.ref} \
                -O BAM \
                -o {output.BAM} \
                {input.SAM}
            
            {EXEC['SAMTOOLS']} index {output.BAM}
            """
        )
#


rule extract_barcodes:
    input:
        BAM = "{OUTDIR}/{SAMPLE}/ont/sorted.bam",
        BB_WHITELIST = "{OUTDIR}/{SAMPLE}/bb/whitelist.txt"
    output:
        BAM = "{OUTDIR}/{SAMPLE}/ont/sorted.bam", #temp()
        TSV = "{OUTDIR}/{SAMPLE}/ont/bc_counts.tsv",
    params:
        KIT = "3prime", #['3prime', '5prime', 'multiome']
        adapter1_suff_length = config["BARCODE_ADAPTER1_SUFF_LENGTH"],
        barcode_min_qv = config["BARCODE_MIN_QUALITY"],
        barcode_length = lambda w: get_barcode_length(w),
        umi_length = lambda w: get_umi_length(w),
    threads: 
        config["CORES"]
    # conda:
    #     "../envs/barcodes.yml"
    run:
        shell(
            f"""
            python scripts/py/extract_barcode.py \
                -t {threads} \
                --kit {params.kit} \
                --adapter1_suff_length {params.adapter1_suff_length} \
                --min_barcode_qv {params.barcode_min_qv} \
                --barcode_length {params.barcode_length} \
                --umi_length {params.umi_length} \
                --output_bam {output.BAM} \
                --output_barcodes {output.TSV} \
                {input.BAM} {input.BB_WHITELIST}
            """
        )
#

# rule cleanup_headers_1:
#     input:
#         BAM_BC_UNCORR_TMP,
#     output:
#         bam=temp(BAM_BC_UNCORR),
#         bai=temp(BAM_BC_UNCORR_BAI),
#     # conda:
#     #     "../envs/{EXEC['SAMTOOLS']}.yml"
#     run:
#         shell(
#             f"""
#             {EXEC['SAMTOOLS']} reheader \
#                 --no-PG \
#                 -c 'grep -v ^@PG' \
#                 {input} \
#             > {output.bam}
            
#             {EXEC['SAMTOOLS']} index {output.bam}
#             """
#         )


# rule assign_barcodes:
#     input:
#         bam=CHROM_BAM_BC_UNCORR,
#         bai=CHROM_BAM_BC_UNCORR_BAI,
#         whitelist=BARCODE_WHITELIST,
#     output:
#         bam=temp(CHROM_BAM_BC_TMP),
#         counts=CHROM_ASSIGNED_BARCODE_COUNTS,
#     params:
#         max_ed=config["BARCODE_MAX_ED"],
#         min_ed_diff=config["BARCODE_MIN_ED_DIFF"],
#         # kit=lambda w: sample_sheet.loc[w.run_id, "kit_name"],
#         KIT = "3prime", #['3prime', '5prime', 'multiome']
#         adapter1_suff_length=config["BARCODE_ADAPTER1_SUFF_LENGTH"],
#         barcode_length=lambda w: get_barcode_length(w),
#         umi_length=lambda w: get_umi_length(w),
#         # barcode_length=config["READ_STRUCTURE_BARCODE_LENGTH"],
#         # umi_length=config["READ_STRUCTURE_UMI_LENGTH"],
#     threads: 1
#     # conda:
#     #     "../envs/barcodes.yml"
#     run:
#         shell(
#             f"""
#             touch {input.bai}

#             python scripts/py/assign_barcodes.py \
#                 -t {threads} \
#                 --output_bam {output.bam} \
#                 --output_counts {output.counts} \
#                 --max_ed {params.max_ed} \
#                 --min_ed_diff {params.min_ed_diff} \
#                 --kit {params.kit} \
#                 --adapter1_suff_length {params.adapter1_suff_length} \
#                 --barcode_length {params.barcode_length} \
#                 --umi_length {params.umi_length} \
#                 {input.bam} {input.whitelist}
#             """
#         )


# Generate count matrix w/ umi-tools for miRNAs
# rule umitools_count_ONT:
#     input:
#         BAM = '{OUTDIR}/{SAMPLE}/ont/aligned.sorted.tagged.bam',
#         BAI = '{OUTDIR}/{SAMPLE}/ont/aligned.sorted.tagged.bam.bai'
#     output:        
#         COUNTS = '{OUTDIR}/{SAMPLE}/ont/counts.tsv.gz'
#     params:
#         OUTDIR = config['OUTDIR']
#     log:
#         log = '{OUTDIR}/{SAMPLE}/ont/count.log'
#     threads:
#         1
#     run:
#         shell(
#             f"""
#             {EXEC['UMITOOLS']} count \
#                 --extract-umi-method=tag \
#                 --per-gene \
#                 --per-cell \
#                 --cell-tag=CB \
#                 --gene-tag=GN \
#                 --umi-tag=UB \
#                 --log={log.log} \
#                 -I {input.BAM} \
#                 -S {output.COUNTS}
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