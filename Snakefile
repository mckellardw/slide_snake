# slide_snake 
## Snakemake workflow to align and quantify spatial transriptomics datasets

import pandas as pd
# import scipy.io
# import scipy.sparse


### Config #############################################################################
configfile:"config/config.yaml"

RECIPE_SHEET = pd.read_csv(
    config["RECIPE_SHEET"], 
    na_filter=False,
    index_col=0
) 


### Directory locations ################################################################
TMPDIR = config["TMPDIR"]
OUTDIR = config["OUTDIR"]


### Executables ########################################################################
EXEC = config["EXEC"]

### Variables and references ###########################################################
SAMPLE_SHEET = pd.read_csv(
    config["SAMPLE_SHEET_PATH"], 
    na_filter=False
)
SAMPLE_SHEET.index = SAMPLE_SHEET["sampleID"]

SAMPLES = list(SAMPLE_SHEET["sampleID"])

# short-read data
R1_FQS = dict(zip(SAMPLES, list(SAMPLE_SHEET["fastq_R1"])))
R1_FQS = {SAMP: READ.split() for SAMP, READ in R1_FQS.items() if READ}
R2_FQS = dict(zip(SAMPLES, list(SAMPLE_SHEET["fastq_R2"])))
R2_FQS = {SAMP: READ.split() for SAMP, READ in R2_FQS.items() if READ}

# long-read data
ONT = dict(zip(SAMPLES, list(SAMPLE_SHEET["ONT"])))
ONT = {SAMP: READ.split() for SAMP, READ in ONT.items() if READ}


### Pre-run setup ######################################################################
# Build dictionaries of recipes & species to use for alignment
RECIPE_DICT = {}        # Dictionary of lists; recipes to use for each sample (short_read module)
RECIPE_ONT_DICT = {}    # Dictionary of lists; recipes to use for each sample (ONT module)
for i in range(0,SAMPLE_SHEET.shape[0]):
    tmp_sample = list(SAMPLE_SHEET["sampleID"])[i]

    # short-read-specific dicts
    if tmp_sample in R2_FQS.keys():
        RECIPE_DICT[tmp_sample] = list(SAMPLE_SHEET["recipe"])[i].split()
    
    # ONT-specific dicts
    if tmp_sample in ONT.keys():
        if len(ONT[tmp_sample]) > 0:
            RECIPE_ONT_DICT[tmp_sample] = list(SAMPLE_SHEET["recipe_ONT"])[i].split()


### recipe_sheet & sample_sheet checks ######################################################################
include: "rules/0_utils.smk" 
check_recipe_sheet(RECIPE_SHEET, RECIPE_DICT, RECIPE_ONT_DICT)
#TODO check_sample_sheet(SAMPLE_SHEET)

### Wildcard constraints ###############################################################
wildcard_constraints:
    OUTDIR = config["OUTDIR"],
    SAMPLE = "[A-Za-z0-9_-]+"


### include rules #######################################################################
# Pre-flight module ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
include: "rules/0a_barcode_maps.smk"

# Short-read module ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## fastq preprocessing & QC
include: "rules/short_read/1a_mergefqs.smk"
include: "rules/short_read/1b_trimming.smk"
include: "rules/short_read/1c_fastqc.smk"

## rRNA Filtering 
include: "rules/short_read/2a_rRNA_bwa.smk"
include: "rules/short_read/2b_rRNA_STAR.smk"
include: "rules/short_read/2c_rRNA_qualimap.smk"
include: "rules/short_read/2d_ribodetector.smk"

## STAR alignment, QC, and post-processing - TODO update numbering
include: "rules/short_read/3a_star_align.smk"
include: "rules/short_read/3b_star_unmapped.smk"
include: "rules/short_read/3c_star_dedup.smk"
include: "rules/short_read/3d_star_qualimap.smk"

## kallisto/bustools alignment
# include: "rules/short_read/4a_kallisto.smk"
include: "rules/short_read/4a_kbpython.smk"
# include: "rules/short_read/4b_kallisto_pseudobam.smk"
# include: "rules/short_read/4c_kallisto_velo.smk"

## small RNA stuff
# include: "rules/short_read/5_mirge.smk"

# ONT module ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## preprocessing
include: "rules/ont/1a_preprocessing.smk"
include: "rules/ont/1b_trimming.smk"
include: "rules/ont/1c_barcode_calling.smk"

## alignment
include: "rules/ont/1d_minimap2.smk"
include: "rules/ont/1d_STAR.smk"

## QC
include: "rules/ont/2_qualimap.smk"
include: "rules/ont/2_read_qc.smk"
include: "rules/ont/2_fastqc.smk"

# Final outputs module ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## scanpy stuff
include: "rules/6a_scanpy_init.smk"

### target rule(s) #####################################################################
rule all:
    input:
        ### short-read targets ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # Module 1 - trimming & QC
        [f"{OUTDIR}/{SAMPLE}/fastqc/{TRIM}_{READ}"
            for SAMPLE in R2_FQS.keys()
            for TRIM in ["preCutadapt","postCutadapt","twiceCutadapt","rRNA_bwa","rRNA_STAR"] 
            for READ in ["R1","R2"] 
        ],  # fastQC results        
        
        # # Module 2 - rRNA filtering        
        # [f"{OUTDIR}/{SAMPLE}/qualimap/rRNA/{TOOL}/{FILE}"
        #     for SAMPLE in R2_FQS.keys() 
        #     for TOOL in ["bwa","STAR"]
        #     for FILE in ["qualimapReport.html","rnaseq_qc_results.csv"] 
        # ], # alignment QC with qualimap [rRNA alignments]
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
        [f"{OUTDIR}/{SAMPLE}/STARsolo/short_read/{RECIPE}/Aligned.sortedByCoord.dedup.out{STRAND}.bam.bai"
            for SAMPLE in R2_FQS.keys() 
            for RECIPE in RECIPE_DICT[SAMPLE] 
            for STRAND in ["", ".fwd", ".rev"]
        ], # deduped and/or strand-split, umi_tools deduplicated .bam #TODO- REF=["STARsolo_rRNA", "STARsolo"]
        [f"{OUTDIR}/{SAMPLE}/STARsolo/short_read/{RECIPE}/Solo.out/GeneFull/raw/matrix.mtx.gz" 
            for SAMPLE in R2_FQS.keys() 
            for RECIPE in RECIPE_DICT[SAMPLE]
        ], # STAR count mats             
        [f"{OUTDIR}/{SAMPLE}/qualimap/{TOOL}/{RECIPE}/{FILE}"
            for SAMPLE in R2_FQS.keys() 
            for RECIPE in RECIPE_DICT[SAMPLE] 
            for TOOL in ["STAR"]
            for FILE in ["qualimapReport.html","rnaseq_qc_result.csv"] 
        ], # alignment QC with qualimap | requires deduped input!    

        [f"{OUTDIR}/{SAMPLE}/fastqc/unmapped/{RECIPE}" 
            for SAMPLE in R2_FQS.keys() 
            for RECIPE in RECIPE_DICT[SAMPLE]
        ], #fastQC results for unmapped reads
        # [f"{OUTDIR}/{SAMPLE}/unmapped/{RECIPE}/blast/Unmapped.out.mate2_blastResults.txt"
        #     for SAMPLE in R2_FQS.keys()
        #     for RECIPE in RECIPE_DICT[SAMPLE]
        # ], # Top BLAST results for unmapped R2 reads   
        
        # Module 4 - kallisto & bustools
        #TODO rewrite for kbpython
        # expand( # kallisto/bustools count mats #TODO
        #     "{OUTDIR}/{SAMPLE}/kb_velo/{LAYER}/output.mtx.gz",
        #     OUTDIR=config["OUTDIR"],
        #     LAYER=["spliced", "unspliced"],
        #     SAMPLE=R2_FQS.keys()
        # ),

        # Module 5 - small RNA
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
        [f"{OUTDIR}/{SAMPLE}/kbpython_{KB}/{RECIPE}/counts_unfiltered/output.h5ad" 
            for SAMPLE in R2_FQS.keys() 
            for RECIPE in RECIPE_DICT[SAMPLE] 
            for KB in ["std"] #TODO "nac", "tcc" 
        ], # anndata files (with spatial info) - kallisto 

        
        ### ONT targets ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        [f"{OUTDIR}/{SAMPLE}/ont/{FILE}" 
            for SAMPLE in ONT.keys() 
            for RECIPE in RECIPE_ONT_DICT[SAMPLE]
            for FILE in [
                    f"minimap2/{RECIPE}/sorted.bam",
                    f"minimap2/{RECIPE}/sorted_gn_cb.bam",
                    f"barcodes_umis/{RECIPE}/read_barcodes_corrected.tsv", #_corrected
                    f"minimap2/{RECIPE}/raw/output.h5ad",
                ]
        ], # ONT outputs

        [f"{OUTDIR}/{SAMPLE}/fastqc/{TRIM}"
            for SAMPLE in ONT.keys()
            for READ in ["R1","R2"]
            for TRIM in ["ont_preAdapterScan",f"ont_preCutadapt_{READ}",f"ont_postCutadapt_{READ}"]
        ], # ONT fastqc - not really useful...

        [f"{OUTDIR}/{SAMPLE}/ont/readqc/{TRIM}_qc.{FILE}"
            for SAMPLE in ONT.keys()
            for READ in ["R1","R2"]
            for RECIPE in RECIPE_ONT_DICT[SAMPLE]
            for TRIM in [f"0_rawInput/merged", f"1_preCutadapt/{READ}", f"2_postCutadapt/{READ}", f"3_aligned/{RECIPE}"]
            for FILE in ["tsv", "png"]
        ], # ONT readqc

        [f"{OUTDIR}/{SAMPLE}/qualimap/ont/{TOOL}/{RECIPE}/{FILE}"
            for SAMPLE in ONT.keys() 
            for RECIPE in RECIPE_ONT_DICT[SAMPLE]
            for TOOL in ["minimap2",]# f"STARsolo/{RECIPE}"
            for FILE in ["qualimapReport.html","rnaseq_qc_results.csv"] 
        ], # alignment QC with qualimap

        # [f"{OUTDIR}/{SAMPLE}/STARsolo/ont/{RECIPE}/{FILE}" 
        #     for SAMPLE in ONT.keys() 
        #     for RECIPE in RECIPE_DICT[SAMPLE]
        #     for FILE in [
        #         ]
        # ], # ONT outputs

        ### DEPRECATED ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # [f"{OUTDIR}/{SAMPLE}/{KB}/{RECIPE}/raw/output.h5ad" 
        #     for SAMPLE in R2_FQS.keys() 
        #     for RECIPE in RECIPE_DICT[SAMPLE] 
        #     for KB in ["kb"] # "kb_velo", "kb_nuc" 
        # ], 
#end

