# slide_snake 
## Snakemake workflow to align and quantify spatial transriptomics datasets

import pandas as pd
import scipy.io
import scipy.sparse

### Config #############################################################################
configfile:"config/config.yaml"

RECIPE_SHEET = pd.read_csv(
    # config["RECIPE_SHEET"], 
    "resources/recipe_sheet.csv",
    na_filter=False,
    index_col=0
) 

### Directory locations ################################################################
TMPDIR = config["TMPDIR"]
OUTDIR = config["OUTDIR"]

### Conda environment locations ########################################################

### Variables and references ###########################################################
SAMPLE_SHEET = pd.read_csv(
    config["SAMPLE_SHEET_PATH"], 
    na_filter=False
)

SAMPLES = list(SAMPLE_SHEET["sampleID"])

# short-read data
R1_FQS = dict(zip(SAMPLES, list(SAMPLE_SHEET["fastq_R1"])))
R1_FQS = {SAMP: READ.split() for SAMP, READ in R1_FQS.items() if READ}
R2_FQS = dict(zip(SAMPLES, list(SAMPLE_SHEET["fastq_R2"])))
R2_FQS = {SAMP: READ.split() for SAMP, READ in R2_FQS.items() if READ}

# long-read data
ONT = dict(zip(SAMPLES, list(SAMPLE_SHEET["ONT"]))) 
ONT = {SAMP: READ.split() for SAMP, READ in ONT.items() if READ}

### Executables ########################################################################
EXEC = config["EXEC"]

### Pre-run setup ######################################################################
# Build dictionaries of recipes & species to use for alignment
RECIPE_DICT = {}    # Dictionary of recipes to use for each sample
RECIPE_ONT_DICT = {}    # Dictionary of recipes to use for each sample
rRNA_STAR_DICT = {} # Dictionary of rRNA reference genomes to use w/ STAR
rRNA_BWA_DICT = {}  # Dictionary of rRNA reference genomes to use w/ bwa
REF_DICT = {}       # Dictionary of reference genomes to use
GTF_DICT = {}       # Dictionary of gene annotations (.gtf format)
IDX_DICT = {}       # Dictionary of kallisto indices
T2G_DICT = {}       # Dictionary of kallisto transcript-to-gene maps
IDX_VELO_DICT = {}  # Dictionary of kallisto indices for RNA velocity
T2G_VELO_DICT = {}  # Dictionary of kallisto transcript-to-gene maps for RNA velocity
BB_DICT = {}        # Dictionary of bead barcode maps
SPECIES_DICT = {}   # Dictionary of species listed for mirge3 analysis

#TODO- add checks so only needed variables are pulled
for i in range(0,SAMPLE_SHEET.shape[0]):
    tmp_sample = list(SAMPLE_SHEET["sampleID"])[i]

    rRNA_STAR_DICT[tmp_sample] = list(SAMPLE_SHEET["STAR_rRNA_ref"])[i]
    rRNA_BWA_DICT[tmp_sample] = list(SAMPLE_SHEET["bwa_rRNA_ref"])[i]
    REF_DICT[tmp_sample] = list(SAMPLE_SHEET["STAR_ref"])[i]
    GTF_DICT[tmp_sample] = list(SAMPLE_SHEET["genes_gtf"])[i]
    BB_DICT[tmp_sample] = list(SAMPLE_SHEET["BB_map"])[i]
    SPECIES_DICT[tmp_sample] = list(SAMPLE_SHEET["species"])[i]

    # short-read-specific dicts
    if R1_FQS[tmp_sample] is not None:
        RECIPE_DICT[tmp_sample] = list(SAMPLE_SHEET["recipe"])[i].split()
        IDX_DICT[tmp_sample] = list(SAMPLE_SHEET["kb_idx"])[i]
        T2G_DICT[tmp_sample] = list(SAMPLE_SHEET["kb_t2g"])[i]
        # IDX_VELO_DICT[tmp_sample] = list(SAMPLE_SHEET["kb_idx_velo"])[i]
        # T2G_VELO_DICT[tmp_sample] = list(SAMPLE_SHEET["kb_t2g_velo"])[i]

    # ONT-specific dicts
    if ONT[tmp_sample] is not None:
        RECIPE_ONT_DICT[tmp_sample] = list(SAMPLE_SHEET["recipe"])[i].split()


### include rules #######################################################################
# Pre-flight module ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
include: "rules/0a_split_bb.smk" #TODO
include: "rules/0_utils.smk" 

# Short-read module ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## fastq preprocessing & QC
include: "rules/short_read/1a_mergefqs.smk"
include: "rules/short_read/1b_trimQC.smk"
include: "rules/short_read/1d_fq2bam.smk"

## rRNA Filtering 
include: "rules/short_read/2a_rRNA_bwa.smk"
include: "rules/short_read/2b_rRNA_STAR.smk"
include: "rules/short_read/2c_rRNA_qualimap.smk"

## STAR alignment, QC, and post-processing - TODO update numbering
include: "rules/short_read/3a_star_align.smk"
include: "rules/short_read/3b_star_unmapped.smk"
include: "rules/short_read/3c_star_dedup.smk"
include: "rules/short_read/3d_star_qualimap.smk"

## kallisto/bustools alignment
include: "rules/short_read/4a_kallisto.smk"
include: "rules/short_read/4a_kbpython.smk"
include: "rules/short_read/4b_kallisto_pseudobam.smk"
include: "rules/short_read/4c_kallisto_velo.smk"

## small RNA stuff
# include: "rules/short_read/5_mirge.smk"
# include: "rules/short_read/5_piRNA_bowtie2.smk"
# include: "rules/short_read/5_miRNA_bowtie2.smk"

# ONT module ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
include: "rules/ont/0_utils.smk"
include: "rules/ont/1a_preprocessing.smk"
include: "rules/ont/1b_cutadapt.smk"
include: "rules/ont/1c_minimap2.smk"
include: "rules/ont/1d_STAR.smk"
include: "rules/ont/1e_qualimap.smk"
include: "rules/ont/2_fastqc.smk"
include: "rules/ont/2_nanoQC.smk"
# include: "rules/ont/2_nanoplot.smk"

# Final outputs module ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## scanpy stuff
include: "rules/6a_scanpy_init.smk"
# include: "rules/6b_mudata_init.smk"

### target rule(s) #####################################################################
rule all:
    input:
        ### ONT targets ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        [f"{OUTDIR}/{SAMPLE}/ont/{FILE}" 
            for SAMPLE in ONT.keys() 
            for RECIPE in RECIPE_ONT_DICT[SAMPLE]
            for FILE in [
                    # "merged_stranded.fq.gz",
                    f"minimap2/{RECIPE}/sorted.bam",
                    f"minimap2/{RECIPE}/sorted_bc.bam",
                    f"minimap2/{RECIPE}/sorted_bc_gn.bam",
                    # f"minimap2/{RECIPE}/counts.tsv",
                    # f"minimap2/{RECIPE}/raw/umitools_counts.tsv.gz",
                    f"minimap2/{RECIPE}/raw/output.h5ad",
                    # "adapter_scan_readids/full_len_R2.fq.gz"
                ]
        ], # ONT outputs
        # [f"{OUTDIR}/{SAMPLE}/STARsolo/ont/{RECIPE}/{FILE}" 
        #     for SAMPLE in ONT.keys() 
        #     for RECIPE in RECIPE_DICT[SAMPLE]
        #     for FILE in [
        #             f"Aligned.sortedByCoord.out.bam"
        #         ]
        # ], # ONT outputs
        [f"{OUTDIR}/{SAMPLE}/fastqc/{TRIM}"
            for SAMPLE in ONT.keys()
            for READ in ["R1","R2"]
            for TRIM in ["ont_preAdapterScan",f"ont_postCutadapt_{READ}"]
        ], # ONT nanoplot
        [f"{OUTDIR}/{SAMPLE}/nanoqc/{TRIM}/{FILE}" 
            for SAMPLE in ONT.keys()
            for READ in ["R1","R2"]
            for TRIM in ["ont_preAdapterScan",f"ont_postCutadapt_{READ}"]
            for FILE in ["nanoQC.html"]
        ], # read QC with nanoQC

        # [f"{OUTDIR}/{SAMPLE}/nanoplot/{TRIM}/summary.txt" 
        #     for SAMPLE in ONT.keys() 
        #     for RECIPE in RECIPE_ONT_DICT[SAMPLE]
        #     for READ in ["R1","R2"]
        #     for TRIM in ["ont_preAdapterScan",f"ont_postCutadapt_{READ}"]
        # ], # ONT NanoPlot
        [f"{OUTDIR}/{SAMPLE}/qualimap/ont/{TOOL}/{RECIPE}/{FILE}"
            for SAMPLE in ONT.keys() 
            for RECIPE in RECIPE_ONT_DICT[SAMPLE]
            for TOOL in ["minimap2",]# f"STARsolo/{RECIPE}"
            for FILE in ["qualimapReport.html","rnaseq_qc_results.csv"] 
        ], # alignment QC with qualimap

        ### short-read targets ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # Module 1 - trimming & QC
        [f"{OUTDIR}/{SAMPLE}/fastqc/{TRIM}_{READ}"
            for SAMPLE in R2_FQS.keys()
            for TRIM in ["preCutadapt","postCutadapt","twiceCutadapt","rRNA_bwa","rRNA_STAR"] 
            for READ in ["R1","R2"] 
        ],  # fastQC results        

        # Module 2 - rRNA filtering        
        [f"{OUTDIR}/{SAMPLE}/qualimap/rRNA/{TOOL}/{FILE}"
            for SAMPLE in R2_FQS.keys() 
            for TOOL in ["bwa"]#,"STARsolo"
            for FILE in ["qualimapReport.html","rnaseq_qc_results.csv"] 
        ], # alignment QC with qualimap [rRNA alignments]
        # expand( #STAR count mats - rRNA
        #     "{OUTDIR}/{SAMPLE}/rRNA/{ALIGNER}/raw/matrix.mtx.gz",
        #     OUTDIR=config["OUTDIR"],
        #     SAMPLE=R2_FQS.keys(),
        #     ALIGNER=[
        #         "STARsolo/Solo.out/GeneFull"
        #         # "bwa" #TODO
        #     ]
        # ),        

        # Module 3 - STAR alignment
        # expand( # deduped and/or strand-split, umi_tools deduplicated .bam #TODO- REF=["STARsolo_rRNA", "STARsolo"]
        #     "{OUTDIR}/{SAMPLE}/STARsolo/short_read/{RECIPE}/Aligned.sortedByCoord.dedup.out{STRAND}.bam.bai",
        #     OUTDIR=config["OUTDIR"],
        #     SAMPLE=R2_FQS.keys(),
        #     STRAND=["", ".fwd", ".rev"]
        # ),
        [f"{OUTDIR}/{SAMPLE}/fastqc/unmapped/{RECIPE}" 
            for SAMPLE in R2_FQS.keys() 
            for RECIPE in RECIPE_DICT[SAMPLE]
        ], #fastQC results for unmapped reads
        # [f"{OUTDIR}/{SAMPLE}/unmapped/{RECIPE}/blast/Unmapped.out.mate2_blastResults.txt",
        #     SAMPLE=R2_FQS.keys()
        #     for RECIPE in RECIPE_DICT[SAMPLE]
        # ], # Top BLAST results for unmapped R2 reads        
        # [f"{OUTDIR}/{SAMPLE}/qualimap/{RECIPE}/{FILE}"
        #     for SAMPLE in R2_FQS.keys() 
        #     for RECIPE in RECIPE_DICT[SAMPLE] 
        #     for FILE in ["qualimapReport.html","rnaseq_qc_result.csv"] 
        # ], # alignment QC with qualimap | requires deduped input!    
        
        # Module 4 - kallisto & bustools
        # [f"{OUTDIR}/{SAMPLE}/STARsolo/short_read/{RECIPE}/Solo.out/GeneFull/raw/matrix.mtx.gz" 
        #     for SAMPLE in R2_FQS.keys() 
        #     for RECIPE in RECIPE_DICT[SAMPLE]
        # ], # STAR count mats

        # Module 5 - small RNA
        # [f"{OUTDIR}/{SAMPLE}/{SMALL}/{RECIPE}/raw/output.h5ad" 
        #     for SAMPLE in R2_FQS.keys() 
        #     for RECIPE in RECIPE_DICT[SAMPLE] 
        #     for SMALL in ["miRNA","piRNA"]
        # ],# anndata files (with spatial info) - small RNA
        # [f"{OUTDIR}/{SAMPLE}/miRge_bulk/{RECIPE}/annotation.report.html" 
        #     for SAMPLE in R2_FQS.keys() 
        #     for RECIPE in RECIPE_DICT[SAMPLE] 
        # ], # miRge3.0 pseudobulk analysis

        # Module 6 - anndata/scanpy
        [f"{OUTDIR}/{SAMPLE}/STARsolo/short_read/{RECIPE}/Solo.out/{SOLO}/raw/{ALGO}.h5ad" 
            for SAMPLE in R2_FQS.keys() 
            for RECIPE in RECIPE_DICT[SAMPLE] 
            for SOLO in ["Gene","GeneFull"]
            for ALGO in ["UniqueAndMult-EM","matrix"]
        ], # anndata files (with spatial info) - STAR        
        [f"{OUTDIR}/{SAMPLE}/{KB}/{RECIPE}/raw/output.h5ad" 
            for SAMPLE in R2_FQS.keys() 
            for RECIPE in RECIPE_DICT[SAMPLE] 
            for KB in ["kb"] # "kb_velo", "kb_nuc" 
        ], # anndata files (with spatial info) - kallisto #TODO- add kb_velo to `KB`        
        [f"{OUTDIR}/{SAMPLE}/{KB}/{RECIPE}/counts_unfiltered/output.h5ad" 
            for SAMPLE in R2_FQS.keys() 
            for RECIPE in RECIPE_DICT[SAMPLE] 
            for KB in ["kbpython"] # "kb_velo", "kb_nuc" 
        ], # anndata files (with spatial info) - kallisto #TODO- add kb_velo to `KB`


## EXTRANEOUS #######################################################################
        # expand( # count matrices for bowtie2 alignment to small RNA reference(s)
        #     "{OUTDIR}/{SAMPLE}/{SMALL_RNA}/{TYPE}",
        #     OUTDIR=config["OUTDIR"],
        #     SAMPLE=R2_FQS.keys(),
        #     SMALL_RNA=["piRNA","miRNA"],
        #     TYPE=["counts.tsv.gz","raw/matrix.mtx.gz"]
        # ),
        # expand( #non-deduplicated .bam
        #     "{OUTDIR}/{SAMPLE}/{REF}/Aligned.sortedByCoord.out.bam.bai",
        #     OUTDIR=config["OUTDIR"],
        #     SAMPLE=R2_FQS.keys(),
        #     REF=["STARsolo_rRNA", "STARsolo"]
        # ),
        # expand( # kallisto/bustools count mats
        #     "{OUTDIR}/{SAMPLE}/kb/raw/output.mtx.gz",
        #     OUTDIR=config["OUTDIR"],
        #     SAMPLE=R2_FQS.keys()
        # ),

        # expand( # kallisto/bustools count mats
        #     "{OUTDIR}/{SAMPLE}/kb_velo/{LAYER}/output.mtx.gz",
        #     OUTDIR=config["OUTDIR"],
        #     LAYER=["spliced", "unspliced"],
        #     SAMPLE=R2_FQS.keys()
        # ),