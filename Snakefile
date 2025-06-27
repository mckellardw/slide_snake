# slide_snake
## Snakemake workflow to align and quantify spatial transriptomics datasets
import pandas as pd
import os


### Config #############################################################################
configfile: "config.yaml"


RECIPE_SHEET = pd.read_csv(config["RECIPE_SHEET"], na_filter=False, index_col=0)

### Directory locations ################################################################
TMPDIR = os.path.abspath(config["TMPDIR"])
OUTDIR = config["OUTDIR"]
# OUTDIR = os.path.abspath(config["OUTDIR"])


### Variables and references ###########################################################
SAMPLE_SHEET = pd.read_csv(config["SAMPLE_SHEET_PATH"], na_filter=False)
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
## Dictionary of lists; recipes to use for each sample (short_read module)
RECIPE_DICT = {}

## Dictionary of lists; recipes to use for each sample (ONT module)
RECIPE_ONT_DICT = {}

for i in range(0, SAMPLE_SHEET.shape[0]):
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
check_sample_sheet(SAMPLE_SHEET)


### Wildcard constraints ###############################################################
wildcard_constraints:
    OUTDIR=config["OUTDIR"],
    SAMPLE="[A-Za-z0-9_-]+",
    RECIPE="[A-Za-z0-9_-]+",


### include rules #######################################################################
# Barcode handling ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
include: "rules/0a_barcode_maps.smk"
# Short-read module ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## fastq preprocessing & QC
include: "rules/short_read/1a_mergefqs.smk"
include: "rules/short_read/1b_trimming.smk"
include: "rules/short_read/1c_barcode_calling.smk"
## rRNA Filtering
include: "rules/short_read/2a_rRNA_bwa.smk"
include: "rules/short_read/2b_ribodetector.smk"
include: "rules/short_read/2c_rRNA_qualimap.smk"
## STAR alignment, QC, and post-processing - TODO update numbering
include: "rules/short_read/3a_star_align.smk"
include: "rules/short_read/3b_star_unmapped.smk"
include: "rules/short_read/3c_star_dedup.smk"
include: "rules/short_read/3q_star_qualimap.smk"
include: "rules/short_read/3u_star_uTAR.smk"
## kallisto/bustools alignment
include: "rules/short_read/4a_kbpython.smk"
## small RNA stuff #TODO
include: "rules/short_read/5a_mirge.smk"
## scanpy stuff
include: "rules/short_read/6a_scanpy_init.smk"
include: "rules/short_read/6b_seurat_init.smk"
include: "rules/short_read/7a_fastqc.smk"
include: "rules/short_read/7b_readqc.smk"
# ONT module ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## preprocessing
include: "rules/ont/1a_preprocessing.smk"
include: "rules/ont/1b_trimming.smk"
include: "rules/ont/1c_barcode_calling.smk"
## alignment
include: "rules/ont/1d_minimap2_genome.smk"
include: "rules/ont/1d_minimap2_transcriptome.smk"
# include: "rules/ont/1e_kallisto-lr.smk"
include: "rules/ont/1f_ultra_genome.smk"
include: "rules/ont/1g_isoquant.smk"
## QC
include: "rules/ont/2a_readqc.smk"
include: "rules/ont/2b_qualimap.smk"


### Build targets #################################################################################

### short-read targets ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Module 1 - preprocessing
## fastQC results
ilmn_fastqc = [
    f"{OUTDIR}/{SAMPLE}/short_read/fastqc/{TRIM}_{READ}"
    for SAMPLE in R2_FQS.keys()
    for TRIM in ["preCutadapt", "postCutadapt", "twiceCutadapt"]  # ,"rRNA_bwa","rRNA_STAR"
    for READ in ["R1", "R2"]
]

ilmn_barcodes = [
    f"{OUTDIR}/{SAMPLE}/short_read/{FILE}"
    for SAMPLE in R2_FQS.keys()
    for RECIPE in RECIPE_DICT[SAMPLE]
    for FILE in [
        f"barcodes_umis/{RECIPE}/barcodes_corrected.tsv",
        f"barcodes_umis/{RECIPE}/bc_correction_stats.txt",
    ]
]

# Module 2 - rRNA filtering
## alignment QC with qualimap [rRNA alignments]
ilmn_rRNA_qualimap = [
    f"{OUTDIR}/{SAMPLE}/short_read/qualimap/rRNA/{TOOL}/{FILE}"
    for SAMPLE in R2_FQS.keys()
    for TOOL in ["bwa"]
    for FILE in [
        "rnaseq/report.pdf",
        "rnaseq/rnaseq_qc_results.csv",
        "bamqc/report.pdf",
        "bamqc/genome_results.csv",
    ]
]


# Module 3 - STAR alignment
## deduped and/or strand-split, umi_tools deduplicated .bam

### STAR count mats
ilmn_STAR_counts = [
    f"{OUTDIR}/{SAMPLE}/short_read/STARsolo/{RECIPE}/Solo.out/{SOLO}/raw/{ALGO}.mtx.gz"
    for SAMPLE in R2_FQS.keys()
    for RECIPE in RECIPE_DICT[SAMPLE]
    for SOLO in ["Gene", "GeneFull"]
    for ALGO in ["UniqueAndMult-EM", "matrix"]
]

### Deduplicated and strand-split alignment files
ilmn_STAR_dedup_bams = [
    f"{OUTDIR}/{SAMPLE}/short_read/{REF}/{RECIPE}/Aligned.sortedByCoord.out{STRAND}{STEP}.{FILE}"
    for SAMPLE in R2_FQS.keys()
    for REF in ["STARsolo"]
    for RECIPE in RECIPE_DICT[SAMPLE]
    for STRAND in ["", ".fwd", ".rev"]
    for STEP in ["", ".dedup"]
    for FILE in ["bam", "bam.bai"]  # TODO add bigWigs
]

## alignment QC with qualimap | requires deduped input!
ilmn_STAR_qualimap = [
    f"{OUTDIR}/{SAMPLE}/short_read/qualimap/{TOOL}/{RECIPE}/{DEDUP}/{FILE}"
    for SAMPLE in R2_FQS.keys()
    for TOOL in ["STAR"]
    for RECIPE in RECIPE_DICT[SAMPLE]
    for DEDUP in ["raw"]  # , "dedup"
    for FILE in [
        "rnaseq/report.pdf",
        "rnaseq/rnaseq_qc_results.csv",
        "bamqc/report.pdf",
        "bamqc/genome_results.csv",
    ]
]

## fastQC results for unmapped reads
ilmn_STAR_unmapped_fastqc = [
    f"{OUTDIR}/{SAMPLE}/short_read/fastqc/unmapped/{RECIPE}"
    for SAMPLE in R2_FQS.keys()
    for RECIPE in RECIPE_DICT[SAMPLE]
]

### STAR uTAR outputs
ilmn_STAR_uTAR = [
    f"{OUTDIR}/{SAMPLE}/short_read/STARsolo/{RECIPE}/TAR/{FILE}"
    for SAMPLE in R2_FQS.keys()
    for RECIPE in RECIPE_DICT[SAMPLE]
    for FILE in ["uTAR.mtx.gz", "qc_plots.png"]
]

# Module 4 - kallisto & bustools
##
# TODO

# Module 5 - small RNA
## miRge3.0 pseudobulk analysis
ilmn_mirge_bulk = [
    f"{OUTDIR}/{SAMPLE}/short_read/miRge_bulk/{RECIPE}/annotation.report.html"
    for SAMPLE in R2_FQS.keys()
    for RECIPE in RECIPE_DICT[SAMPLE]
]

# Module 6 - anndata/scanpy
## anndata and seurat/rds files (with spatial info) - STAR
ilmn_STAR_cache = [
    f"{OUTDIR}/{SAMPLE}/short_read/STARsolo/{RECIPE}/Solo.out/{SOLO}/raw/{ALGO}.{FILE}"
    for SAMPLE in R2_FQS.keys()
    for RECIPE in RECIPE_DICT[SAMPLE]
    for SOLO in ["Gene", "GeneFull"]
    for ALGO in ["UniqueAndMult-EM", "matrix"]
    for FILE in ["h5ad", "rds"]
]

## anndata files (with spatial info) - kallisto
ilmn_kb_cache = [
    f"{OUTDIR}/{SAMPLE}/short_read/kbpython_{KB}/{RECIPE}/counts_unfiltered/output.{FILE}"
    for SAMPLE in R2_FQS.keys()
    for RECIPE in RECIPE_DICT[SAMPLE]
    for KB in ["std"]  # TODO "nac", "tcc"
    for FILE in [
        "h5ad",
        # "rds"
    ]
]

# Module 7 - final QC

# ILMN readqc - custom QC scripts
ilmn_readqc = [
    f"{OUTDIR}/{SAMPLE}/short_read/readqc/{TRIM}_qc.{FILE}"
    for SAMPLE in R2_FQS.keys()
    for READ in ["R1", "R2"]
    for RECIPE in RECIPE_DICT[SAMPLE]
    for TRIM in [
        # f"0_rawInput/{READ}",
        # f"1_preCutadapt/{READ}",
        # f"2_postCutadapt/{READ}",
        # f"3_twiceCutadapt/{READ}",
        f"4_aligned/{RECIPE}",
    ]
    for FILE in ["tsv.gz", "png"]
]


### ONT targets ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
ont_barcodes = [
    f"{OUTDIR}/{SAMPLE}/ont/{FILE}"
    for SAMPLE in ONT.keys()
    for RECIPE in RECIPE_ONT_DICT[SAMPLE]
    for FILE in [
        f"barcodes_umis/{RECIPE}/barcodes_corrected.tsv",
        f"barcodes_umis/{RECIPE}/bc_correction_stats.txt",
    ]
]

ont_preprocessing_summary = [
    f"{OUTDIR}/{SAMPLE}/ont/{FILE}"
    for SAMPLE in ONT.keys()
    for RECIPE in RECIPE_ONT_DICT[SAMPLE]
    for FILE in [
        f"logs/1a_adapter_scan_summary.csv",
        f"plots/1a_adapter_scan_summary.png",
        # f"plots/1b_cutadapt_summary.png",
    ]
]

ont_minimap_genome = [
    f"{OUTDIR}/{SAMPLE}/ont/minimap2/{RECIPE}/{FILE}"
    for SAMPLE in ONT.keys()
    for RECIPE in RECIPE_ONT_DICT[SAMPLE]
    for FILE in [
        f"sorted_filtered_cb_ub_gn.bam",
        # "sorted_filtered_cb_ub_gn_pos.bam",
        f"raw/output.h5ad",
        # f"raw/output.rds",
    ]
]

ont_minimap_txome = [
    f"{OUTDIR}/{SAMPLE}/ont/minimap2_txome/{RECIPE}/{FILE}"
    for SAMPLE in ONT.keys()
    for RECIPE in RECIPE_ONT_DICT[SAMPLE]
    for FILE in [
        f"aligned_filtered_sorted_cb_ub.bam",  # input to oarfish
        f"oarfish/P.meta_info.json",
        # f"aligned_gn_cb.bam",
        # f"raw/output.h5ad",
        # f"raw/output.rds",
    ]
]

ont_ultra_genome = [
    f"{OUTDIR}/{SAMPLE}/ont/ultra/{RECIPE}/{FILE}"
    for SAMPLE in ONT.keys()
    for RECIPE in RECIPE_ONT_DICT[SAMPLE]
    for FILE in [
        f"sorted_filtered_cb_ub.bam",
        f"raw/umitools_counts.tsv.gz",
        f"raw/output.h5ad",
        # f"raw/output.rds",
    ]
]

ont_isoquant = [
    f"{OUTDIR}/{SAMPLE}/ont/{ALIGN}/{RECIPE}/{FILE}"
    for SAMPLE in ONT.keys()
    for ALIGN in ["minimap2"]  # , "ultra"
    for RECIPE in RECIPE_ONT_DICT[SAMPLE]
    for FILE in [
        f"sorted_filtered_cb_ub_gn_ig.bam",
        f"isoquant/output.h5ad",
        # f"raw/output.rds",
    ]
]

# ONT readqc - custom QC scripts
ont_readqc = [
    f"{OUTDIR}/{SAMPLE}/ont/readqc/{TRIM}_qc.{FILE}"
    for SAMPLE in ONT.keys()
    for READ in ["R1", "R2"]
    for RECIPE in RECIPE_ONT_DICT[SAMPLE]
    for ALIGN in ["minimap2"]  # , "ultra"
    for TRIM in [
        # f"0_rawInput/merged",
        f"1_preCutadapt/{READ}",
        f"2_postCutadapt/{READ}",
        f"3_aligned/{ALIGN}/{RECIPE}",
    ]
    for FILE in ["tsv.gz", "png"]
]

# alignment QC with qualimap
ont_qualimap = [
    f"{OUTDIR}/{SAMPLE}/ont/qualimap/{TOOL}/{RECIPE}/{FILE}"
    for SAMPLE in ONT.keys()
    for RECIPE in RECIPE_ONT_DICT[SAMPLE]
    for TOOL in [
        "minimap2",
    ]
    for FILE in [
        "rnaseq/report.pdf",
        "rnaseq/rnaseq_qc_results.csv",
        "bamqc/report.pdf",
        "bamqc/genome_results.csv",
    ]
]


# kallisto-lr outputs
# ont_kb = [f"{OUTDIR}/{SAMPLE}/ont/kb/{RECIPE}/{FILE}"
#     for SAMPLE in ONT.keys()
#     for RECIPE in RECIPE_ONT_DICT[SAMPLE]
#     for FILE in [ ]
# ],


### Target rule #################################################################################
rule all:
    input:
        # ilmn_barcodes,
        # ilmn_rRNA_qualimap,
        #ilmn_STAR_dedup_bams,
        ilmn_STAR_counts,
        #ilmn_STAR_qualimap,
        # ilmn_STAR_unmapped_fastqc,
        # ilmn_STAR_uTAR,
        # ilmn_mirge_bulk,  #
        ilmn_STAR_cache,
        #ilmn_kb_cache,
        ilmn_fastqc,
        ilmn_readqc,
        ont_barcodes,
        ont_preprocessing_summary,
        ont_minimap_genome,
        # ont_minimap_txome,
        ont_ultra_genome,
        ont_isoquant,
        ont_readqc,
        #ont_qualimap,
