########################################################################################################
# slide_snake
#   Snakemake workflow to align and quantify spatial transriptomics datasets
########################################################################################################

import pandas as pd
import scipy.io
import scipy.sparse

########################################################################################################
# Config file
########################################################################################################
configfile:'config.yaml'
RECIPE_SHEET = pd.read_csv(
    config["RECIPE_SHEET"], 
    na_filter=False,
    index_col=0
) #"resources/recipe_sheet.csv"

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
EXEC = config['EXEC']

########################################################################################################
# Pre-run setup
########################################################################################################
# Build dictionaries of recipes & species to use for alignment
RECIPE_DICT = {}    # Dictionary of recipes to use for each sample
rRNA_DICT = {}      # Dictionary of rRNA reference genomes to use
REF_DICT = {}       # Dictionary of reference genomes to use
GTF_DICT = {}       # Dictionary of gene annotations (.gtf format)
IDX_DICT = {}       # Dictionary of kallisto indices
T2G_DICT = {}       # Dictionary of kallisto transcript-to-gene maps
IDX_VELO_DICT = {}  # Dictionary of kallisto indices for RNA velocity
T2G_VELO_DICT = {}  # Dictionary of kallisto transcript-to-gene maps for RNA velocity
BB_DICT = {}        # Dictionary of bead barcode maps
SPECIES_DICT = {}   # Dictionary of species listed for mirge3 analysis

for i in range(0,SAMPLE_SHEET.shape[0]):
    tmp_sample = list(SAMPLE_SHEET["sampleID"])[i]
    RECIPE_DICT[tmp_sample] = list(SAMPLE_SHEET["recipe"])[i].split()
    rRNA_DICT[tmp_sample] = list(SAMPLE_SHEET["STAR_rRNA_ref"])[i]
    REF_DICT[tmp_sample] = list(SAMPLE_SHEET["STAR_ref"])[i]
    GTF_DICT[tmp_sample] = list(SAMPLE_SHEET["genes_gtf"])[i]
    IDX_DICT[tmp_sample] = list(SAMPLE_SHEET["kb_idx"])[i]
    T2G_DICT[tmp_sample] = list(SAMPLE_SHEET["kb_t2g"])[i]
    # IDX_VELO_DICT[tmp_sample] = list(SAMPLE_SHEET["kb_idx_velo"])[i]
    # T2G_VELO_DICT[tmp_sample] = list(SAMPLE_SHEET["kb_t2g_velo"])[i]
    BB_DICT[tmp_sample] = list(SAMPLE_SHEET["BB_map"])[i]
    SPECIES_DICT[tmp_sample] = list(SAMPLE_SHEET["species"])[i]

########################################################################################################
rule all:
    input:
        # expand( # count matrices for bowtie2 alignment to small RNA reference(s)
        #     '{OUTDIR}/{SAMPLE}/{SMALL_RNA}/{TYPE}',
        #     OUTDIR=config['OUTDIR'],
        #     SAMPLE=SAMPLES,
        #     SMALL_RNA=['piRNA','miRNA'],
        #     TYPE=["counts.tsv.gz","raw/matrix.mtx.gz"]
        # ),
        # expand( # miRge3.0 pseudobulk analysis
        #     '{OUTDIR}/{SAMPLE}/miRge_bulk/annotation.report.html',
        #     OUTDIR=config['OUTDIR'],
        #     SAMPLE=SAMPLES
        # ),
        [f"{OUTDIR}/{SAMPLE}/STARsolo/{RECIPE}/Solo.out/{SOLO}/raw/UniqueAndMultEM.h5ad" 
            for SAMPLE in SAMPLES 
            for RECIPE in RECIPE_DICT[SAMPLE] 
            for SOLO in ['Gene','GeneFull']
        ], # anndata files (with spatial info) - STAR
        [f"{OUTDIR}/{SAMPLE}/{KB}/{RECIPE}/raw/output.h5ad" 
            for SAMPLE in SAMPLES 
            for RECIPE in RECIPE_DICT[SAMPLE] 
            for KB in ['kb']
        ], # anndata files (with spatial info) - kallisto #TODO- add kb_velo to `KB`
        # [f"{OUTDIR}/{SAMPLE}/{SMALL}/{RECIPE}/raw/output.h5ad" 
        #     for SAMPLE in SAMPLES 
        #     for RECIPE in RECIPE_DICT[SAMPLE] 
        #     for SMALL in ['miRNA','piRNA']
        # ],# anndata files (with spatial info) - small RNA            
        [f"{OUTDIR}/{SAMPLE}/STARsolo/{RECIPE}/Solo.out/GeneFull/raw/matrix.mtx.gz" 
            for SAMPLE in SAMPLES 
            for RECIPE in RECIPE_DICT[SAMPLE]
        ], # STAR count mats
        expand( #STAR count mats - rRNA
            '{OUTDIR}/{SAMPLE}/rRNA/{ALIGNER}/raw/matrix.mtx.gz',
            OUTDIR=config['OUTDIR'],
            SAMPLE=SAMPLES,
            ALIGNER=[
                "STARsolo/Solo.out/GeneFull"
                # "bwa" #TODO
            ]
        ),
        # expand( #non-deduplicated .bam
        #     '{OUTDIR}/{SAMPLE}/{REF}/Aligned.sortedByCoord.out.bam.bai',
        #     OUTDIR=config['OUTDIR'],
        #     SAMPLE=SAMPLES,
        #     REF=["STARsolo_rRNA", "STARsolo"]
        # ),
        # expand( # kallisto/bustools count mats
        #     '{OUTDIR}/{SAMPLE}/kb/raw/output.mtx.gz',
        #     OUTDIR=config['OUTDIR'],
        #     SAMPLE=SAMPLES
        # ),
        # expand( # kallisto/bustools count mats
        #     '{OUTDIR}/{SAMPLE}/kb_velo/{LAYER}/output.mtx.gz',
        #     OUTDIR=config['OUTDIR'],
        #     LAYER=['spliced','unspliced'],
        #     SAMPLE=SAMPLES
        # ),
        # expand(  # alignment QC with qualimap | requires deduped input!
        #     '{OUTDIR}/{SAMPLE}/qualimap/{FILE}',
        #     OUTDIR=config['OUTDIR'],
        #     SAMPLE=SAMPLES,
        #     FILE=["qualimapReport.html","rnaseq_qc_result.csv"]
        # ),
        # expand( # deduped and/or strand-split, umi_tools deduplicated .bam #TODO- REF=["STARsolo_rRNA", "STARsolo"])
        #     '{OUTDIR}/{SAMPLE}/STARsolo/{RECIPE}/Aligned.sortedByCoord.dedup.out{STRAND}.bam.bai',
        #     OUTDIR=config['OUTDIR'],
        #     SAMPLE=SAMPLES,
        #     STRAND=["", ".fwd", ".rev"]
        # ),
        # expand( #fastQC results for unmapped reads
        #     '{OUTDIR}/{SAMPLE}/fastqc/unmapped',
        #     OUTDIR=config['OUTDIR'],
        #     SAMPLE=SAMPLES
        # ),
        # expand( # blastn results for unmapped R2 reads
        #     '{OUTDIR}/{SAMPLE}/Unmapped.out.mate2_blastResults.txt',
        #     OUTDIR=config['OUTDIR'],
        #     SAMPLE=SAMPLES
        # ),
        # expand( # fastQC results
        #     '{OUTDIR}/{SAMPLE}/fastqc/{TRIM}_{READ}',
        #     OUTDIR=config['OUTDIR'],
        #     SAMPLE=SAMPLES,
        #     READ=["R1","R2"],
        #     TRIM = ["preTrim","postTrim"]
        # )
        

# fastq preprocessing & QC
include: "rules/1a_mergefqs.smk"
include: "rules/1b_trimQC.smk"
include: "rules/1c_split_bb.smk"

# STAR alignment, QC, and post-processing
# include: "rules/2a_rRNA_bwa.smk"
include: "rules/2b_rRNA_STAR.smk"
include: "rules/2b_star_align.smk"
include: "rules/2c_star_unmapped.smk"
include: "rules/2d_star_qualimap.smk"
include: "rules/2e_star_dedup.smk"

# kallisto/bustools alignment
include: "rules/3a_kallisto.smk"
include: "rules/3b_kallisto_pseudobam.smk"
include: "rules/3c_kallisto_velo.smk"

# small RNA stuff
# include: "rules/4_mirge.smk"
# include: "rules/4_piRNA_bowtie2.smk"
# include: "rules/4_miRNA_bowtie2.smk"

# scanpy stuff
include: "rules/5a_scanpy_init.smk"

# Post-processing, prep for downstream analyses
#TODO:
# - Add proper saturation estimation (not just saturation based on reads aligning to known genes)
