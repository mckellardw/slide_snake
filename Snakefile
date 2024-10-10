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
include: "rules/short_read/2b_ribodetector.smk"
include: "rules/short_read/2c_rRNA_qualimap.smk"

## STAR alignment, QC, and post-processing - TODO update numbering
include: "rules/short_read/3a_star_align.smk"
include: "rules/short_read/3b_star_unmapped.smk"
include: "rules/short_read/3c_star_dedup.smk"
include: "rules/short_read/3d_star_qualimap.smk"

## kallisto/bustools alignment
include: "rules/short_read/4a_kbpython.smk"
# include: "rules/short_read/4b_kbpython_velo.smk"

## small RNA stuff #TODO
# include: "rules/short_read/5_mirge.smk"

## scanpy stuff
include: "rules/short_read/6a_scanpy_init.smk"

# ONT module ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## preprocessing
include: "rules/ont/1a_preprocessing.smk"
include: "rules/ont/1b_trimming.smk"
include: "rules/ont/1c_barcode_calling.smk"

## alignment
include: "rules/ont/1d_minimap2_genome.smk"
include: "rules/ont/1d_minimap2_transcriptome.smk"
# include: "rules/ont/1d_STAR.smk"
# include: "rules/ont/1e_kallisto-lr.smk"

## QC
include: "rules/ont/2_qualimap.smk"
include: "rules/ont/2_read_qc.smk"
include: "rules/ont/2_fastqc.smk"



### Build targets #################################################################################
### short-read targets ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Module 1 - trimming & QC
ilmn_fastqc = [f"{OUTDIR}/{SAMPLE}/fastqc/{TRIM}_{READ}"
    for SAMPLE in R2_FQS.keys()
    for TRIM in ["preCutadapt","postCutadapt","twiceCutadapt"] #,"rRNA_bwa","rRNA_STAR"
    for READ in ["R1","R2"] 
]  # fastQC results        

# Module 2 - rRNA filtering        
ilmn_rRNA_qualimap = [f"{OUTDIR}/{SAMPLE}/qualimap/rRNA/{TOOL}/{FILE}"
    for SAMPLE in R2_FQS.keys() 
    for TOOL in ["bwa"] #,"STAR"
    for FILE in ["qualimapReport.html","rnaseq_qc_results.csv"] 
] # alignment QC with qualimap [rRNA alignments]

# Module 3 - STAR alignment        
ilmn_STAR_dedup_bams = [f"{OUTDIR}/{SAMPLE}/{REF}/short_read/{RECIPE}/Aligned.sortedByCoord.{STEP}out{STRAND}.bam.bai"
        for SAMPLE in R2_FQS.keys() 
        for REF in ["STARsolo"]
        for RECIPE in RECIPE_DICT[SAMPLE] 
        for STEP in ["", "dedup."]
        for STRAND in ["", ".fwd", ".rev"]
    ] # deduped and/or strand-split, umi_tools deduplicated .bam 

ilmn_STAR_counts = [f"{OUTDIR}/{SAMPLE}/short_read/STARsolo/{RECIPE}/Solo.out/{SOLO}/raw/{ALGO}.mtx.gz" 
    for SAMPLE in R2_FQS.keys() 
    for RECIPE in RECIPE_DICT[SAMPLE]
    for SOLO in ["Gene","GeneFull"]
    for ALGO in ["UniqueAndMult-EM","matrix"]
] # STAR count mats 

ilmn_STAR_dedup_qualimap = [f"{OUTDIR}/{SAMPLE}/qualimap/{TOOL}/{RECIPE}/{FILE}"
    for SAMPLE in R2_FQS.keys() 
    for RECIPE in RECIPE_DICT[SAMPLE] 
    for TOOL in ["STAR"]
    for FILE in ["qualimapReport.html","rnaseq_qc_result.csv"] 
] # alignment QC with qualimap | requires deduped input! 

ilmn_STAR_unmapped_fastqc = [f"{OUTDIR}/{SAMPLE}/fastqc/unmapped/{RECIPE}" 
    for SAMPLE in R2_FQS.keys() 
    for RECIPE in RECIPE_DICT[SAMPLE]
] #fastQC results for unmapped reads

# Module 4 - kallisto & bustools
##

# Module 5 - small RNA
ilmn_mirge_bulk = [f"{OUTDIR}/{SAMPLE}/miRge_bulk/{RECIPE}/annotation.report.html" 
    for SAMPLE in R2_FQS.keys() 
    for RECIPE in RECIPE_DICT[SAMPLE] 
] # miRge3.0 pseudobulk analysis

# Module 6 - anndata/scanpy
ilmn_STAR_h5ad = [f"{OUTDIR}/{SAMPLE}/short_read/STARsolo/{RECIPE}/Solo.out/{SOLO}/raw/{ALGO}.h5ad" 
    for SAMPLE in R2_FQS.keys() 
    for RECIPE in RECIPE_DICT[SAMPLE] 
    for SOLO in ["Gene","GeneFull"]
    for ALGO in ["UniqueAndMult-EM","matrix"]
] # anndata files (with spatial info) - STAR 

ilmn_kb_h5ad = [f"{OUTDIR}/{SAMPLE}/kbpython_{KB}/{RECIPE}/counts_unfiltered/output.h5ad" 
    for SAMPLE in R2_FQS.keys() 
    for RECIPE in RECIPE_DICT[SAMPLE] 
    for KB in ["std"] #TODO "nac", "tcc" 
] # anndata files (with spatial info) - kallisto 


### ONT targets ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
ont_minimap = [f"{OUTDIR}/{SAMPLE}/ont/{FILE}" 
    for SAMPLE in ONT.keys() 
    for RECIPE in RECIPE_ONT_DICT[SAMPLE]
    for FILE in [
            f"minimap2/{RECIPE}/sorted_gn_cb.bam",
            f"barcodes_umis/{RECIPE}/read_barcodes_corrected.tsv",
            f"barcodes_umis/{RECIPE}/bc_correction_stats.txt",
            f"minimap2/{RECIPE}/raw/output.h5ad",
            # f"minimap2_txome/{RECIPE}/raw/output.h5ad",
        ]
] # ONT outputs

ont_fastqc = [f"{OUTDIR}/{SAMPLE}/fastqc/{TRIM}"
    for SAMPLE in ONT.keys()
    for READ in ["R1","R2"]
    for TRIM in ["ont_preAdapterScan", f"ont_preCutadapt_{READ}", f"ont_postCutadapt_{READ}"]
] # ONT fastqc - not really useful, but I coded it out...

ont_readqc = [f"{OUTDIR}/{SAMPLE}/ont/readqc/{TRIM}_qc.{FILE}"
    for SAMPLE in ONT.keys()
    for READ in ["R1","R2"]
    for RECIPE in RECIPE_ONT_DICT[SAMPLE]
    for TRIM in [f"0_rawInput/merged", f"1_preCutadapt/{READ}", f"2_postCutadapt/{READ}", f"3_aligned/{RECIPE}"]
    for FILE in ["tsv", "png"]
] # ONT readqc

ont_qualimap = [f"{OUTDIR}/{SAMPLE}/qualimap/ont/{TOOL}/{RECIPE}/{FILE}"
    for SAMPLE in ONT.keys() 
    for RECIPE in RECIPE_ONT_DICT[SAMPLE]
    for TOOL in ["minimap2",]# f"STARsolo/{RECIPE}"
    for FILE in ["qualimapReport.html", "rnaseq_qc_results.csv"] 
], # alignment QC with qualimap


#ont_kb = [f"{OUTDIR}/{SAMPLE}/ont/kb/{RECIPE}/{FILE}" 
#     for SAMPLE in ONT.keys() 
#     for RECIPE in RECIPE_DICT[SAMPLE]
#     for FILE in [ ]
# ], # ONT outputs


### Target rule #################################################################################
rule all:
    input:
        ont_minimap,
        # ont_readqc,
        ont_qualimap,
        # ont_fastqc,

        ilmn_fastqc,
        # ilmn_rRNA_qualimap,
        # ilmn_STAR_dedup_bams,
        # ilmn_STAR_counts,
        # ilmn_STAR_dedup_qualimap,
        # ilmn_STAR_unmapped_fastqc,
        # ilmn_mirge_bulk,
        ilmn_STAR_h5ad,
        # ilmn_kb_h5ad,
