########################################################################################################
# slide_snake
#   Snakemake workflow to align and quantify Seeker/SlideSeq datasets
#   Written by David McKellar
########################################################################################################

import pandas as pd
import scipy.io
import scipy.sparse

########################################################################################################
# Config file
########################################################################################################
configfile:'config.yaml'
RECIPE_SHEET = pd.read_csv(config["RECIPE_SHEET"], na_filter=False,index_col=0) #"resources/recipe_sheet.csv"
########################################################################################################
# Directories and locations
########################################################################################################
TMPDIR = config['TMPDIR']
OUTDIR = config['OUTDIR']

########################################################################################################
# Variables and references
########################################################################################################
SAMPLE_SHEET = pd.read_csv(config["SAMPLE_SHEET_PATH"], na_filter=False)
SAMPLES = list(SAMPLE_SHEET['sampleID'])

R1_FQS = dict(zip(SAMPLES, list(SAMPLE_SHEET['fastq_R1'])))
R2_FQS = dict(zip(SAMPLES, list(SAMPLE_SHEET['fastq_R2'])))

########################################################################################################
# Executables
########################################################################################################
STAR_EXEC = config["STAR_EXEC"]
FASTQC_EXEC = config["FASTQC_EXEC"]
CUTADAPT_EXEC = config["CUTADAPT_EXEC"]
SAMTOOLS_EXEC = config["SAMTOOLS_EXEC"]
UMITOOLS_EXEC = config["UMITOOLS_EXEC"]
QUALIMAP_EXEC = config["QUALIMAP_EXEC"]

########################################################################################################
# Pre-run setup
########################################################################################################
# Build dictionaries of recipes & species to use for alignment
RECIPE_DICT = {} # Dictionary of recipes to use for each sample
rRNA_DICT = {} # Dictionary of rRNA reference genomes to use
REF_DICT = {} # Dictionary of reference genomes to use
GTF_DICT = {} # Dictionary of gene annotations (.gtf format)
IDX_DICT = {} # Dictionary of kallisto indices
T2G_DICT = {} # Dictionary of kallisto transcript-to-gene maps
BB_DICT = {} # Dictionary of bead barcode maps
for i in range(0,SAMPLE_SHEET.shape[0]):
    tmp_sample = list(SAMPLE_SHEET["sampleID"])[i]
    RECIPE_DICT[tmp_sample] = list(SAMPLE_SHEET["recipe"])[i] 
    rRNA_DICT[tmp_sample] = list(SAMPLE_SHEET["STAR_rRNA_ref"])[i]
    REF_DICT[tmp_sample] = list(SAMPLE_SHEET["STAR_ref"])[i]
    GTF_DICT[tmp_sample] = list(SAMPLE_SHEET["genes_gtf"])[i]
    IDX_DICT[tmp_sample] = list(SAMPLE_SHEET["kb_idx"])[i]
    T2G_DICT[tmp_sample] = list(SAMPLE_SHEET["kb_t2g"])[i]
    BB_DICT[tmp_sample] = list(SAMPLE_SHEET["BB_map"])[i]

########################################################################################################
rule all:
    input:
        # expand('{OUTDIR}/{sample}/{REF}/Solo.out/GeneFull/raw/UniqueAndMultEM.h5ad', OUTDIR=config['OUTDIR'], sample=SAMPLES, REF=["STARsolo_rRNA", "STARsolo"]), 
        # expand('{OUTDIR}/{sample}/{REF}/Solo.out/GeneFull/raw/barcodes_noUnderscore.tsv.gz', OUTDIR=config['OUTDIR'], sample=SAMPLES, REF=["STARsolo_rRNA", "STARsolo"]), #Barcode lists w/ underscores removed
        expand( #STAR count mats
            '{OUTDIR}/{sample}/{REF}/Solo.out/GeneFull/raw/matrix.mtx.gz', 
            OUTDIR=config['OUTDIR'], 
            sample=SAMPLES, 
            REF=["STARsolo_rRNA", "STARsolo"]
            ), 
        # expand('{OUTDIR}/{sample}/kb/counts_unfiltered/output.mtx.gz', OUTDIR=config['OUTDIR'], sample=SAMPLES), #kallisto count mats
        expand(  # alignment QC qith qualimap | requires deduped input!
            '{OUTDIR}/{sample}/qualimap/qualimapReport.html', 
            OUTDIR=config['OUTDIR'], 
            sample=SAMPLES
            ),
        expand( #TODO- REF=["STARsolo_rRNA", "STARsolo"]), # umi_tools deduplicated .bam
            '{OUTDIR}/{sample}/STARsolo/Aligned.sortedByCoord.dedup.out.bam.bai', 
            OUTDIR=config['OUTDIR'], 
            sample=SAMPLES
        ), 
        expand(
            '{OUTDIR}/{sample}/STARsolo/Aligned.sortedByCoord.dedup.out.{STRAND}.bam.bai', OUTDIR=config['OUTDIR'], sample=SAMPLES, STRAND=["fwd", "rev"]), # umi_tools deduplicated .bam
        # expand('{OUTDIR}/{sample}/{REF}/Aligned.sortedByCoord.dedup.out_plus.bw', OUTDIR=config['OUTDIR'], sample=SAMPLES, REF=["STARsolo_rRNA", "STARsolo"]), # strand-split bigWigs
        expand( # fastQC results
            '{OUTDIR}/{sample}/{TRIM}_fastqc_{READ}', 
            OUTDIR=config['OUTDIR'],
            sample=SAMPLES, 
            READ=["R1","R2"], 
            TRIM = ["preTrim","postTrim"]
        ),  
        expand( #non-deduplicated .bam; used for saturation estimation
            '{OUTDIR}/{sample}/{REF}/Aligned.sortedByCoord.out.bam.bai', 
            OUTDIR=config['OUTDIR'], 
            sample=SAMPLES, 
            REF=["STARsolo_rRNA", "STARsolo"]
        ), 
        expand( #fastQC results for unmapped reads
            '{OUTDIR}/{sample}/Unmapped_fastqc', 
            OUTDIR=config['OUTDIR'], 
            sample=SAMPLES
        ), 
        expand( # blastn results for unmapped R2 reads 
            '{OUTDIR}/{sample}/Unmapped.out.mate2_blastResults.txt', 
            OUTDIR=config['OUTDIR'], 
            sample=SAMPLES
        ), 

# fastq preprocessing & QC
include: "rules/1a_mergefqs.smk"
include: "rules/1b_trimQC.smk"
include: "rules/1c_split_bb.smk"

# STAR alignment, QC, and post-processing
include: "rules/2a_star_align_rRNA.smk"
include: "rules/2b_star_align.smk"
include: "rules/2c_star_unmapped.smk"
include: "rules/2d_star_qualimap.smk"
include: "rules/2e_star_dedup.smk"

# kallisto/bustools alignment
include: "rules/3a_kallisto_align.smk"
include: "rules/3b_kallisto_pseudobam.smk"

# scanpy stuff
include: "rules/4a_scanpy_init.smk"

# Post-processing, prep for downstream analyses
#TODO:
# - Initialize a .h5ad object for easy loading into python later
#   - Add spatial location!

