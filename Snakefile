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
CHEMISTRY_SHEET = pd.read_csv(config["CHEMISTRY_SHEET"], na_filter=False,index_col=0) #"resources/chemistry_sheet.csv"
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
# Build dictionaries of chemistries & species to use for alignment
CHEM_DICT = {} # Dictionary of chemistry recipe to use for each sample
REF_DICT = {} # Dictionary of reference genomes to use
GTF_DICT = {} # Dictionary of gene annotations (.gtf format)
IDX_DICT = {} # Dictionary of kallisto indices
T2G_DICT = {} # Dictionary of kallisto transcript-to-gene maps
BB_DICT = {} # Dictionary of bead barcode maps
for i in range(0,SAMPLE_SHEET.shape[0]):
    tmp_sample = list(SAMPLE_SHEET["sampleID"])[i]
    CHEM_DICT[tmp_sample] = list(SAMPLE_SHEET["chemistry"])[i]
    REF_DICT[tmp_sample] = list(SAMPLE_SHEET["STAR_ref"])[i]
    GTF_DICT[tmp_sample] = list(SAMPLE_SHEET["genes_gtf"])[i]
    IDX_DICT[tmp_sample] = list(SAMPLE_SHEET["kb_idx"])[i]
    T2G_DICT[tmp_sample] = list(SAMPLE_SHEET["kb_t2g"])[i]
    BB_DICT[tmp_sample] = list(SAMPLE_SHEET["BB_map"])[i]

########################################################################################################
rule all:
    input:
        # expand('{OUTDIR}/{sample}/{REF}/Solo.out/Gene/raw/UniqueAndMultEM.h5ad', OUTDIR=config['OUTDIR'], sample=SAMPLES, REF=["STARsolo_rRNA", "STARsolo"]), 
        # expand('{OUTDIR}/{sample}/{REF}/Solo.out/Gene/raw/barcodes_noUnderscore.tsv.gz', OUTDIR=config['OUTDIR'], sample=SAMPLES, REF=["STARsolo_rRNA", "STARsolo"]), #Barcode lists w/ underscores removed
        expand('{OUTDIR}/{sample}/{REF}/Solo.out/Gene/raw/matrix.mtx.gz', OUTDIR=config['OUTDIR'], sample=SAMPLES, REF=["STARsolo_rRNA", "STARsolo"]), #STAR count mats
        expand('{OUTDIR}/{sample}/kb/counts_unfiltered/output.mtx.gz', OUTDIR=config['OUTDIR'], sample=SAMPLES), #kallisto count mats
        expand('{OUTDIR}/{sample}/qualimap/qualimapReport.html', OUTDIR=config['OUTDIR'], sample=SAMPLES), # alignment QC qith qualimap | requires deduped input!
        # expand('{OUTDIR}/{sample}/{REF}/Aligned.sortedByCoord.dedup.out.bam.bai', OUTDIR=config['OUTDIR'], sample=SAMPLES, REF=["STARsolo_rRNA", "STARsolo"]), # umi_tools deduplicated .bam
        # expand('{OUTDIR}/{sample}/{REF}/Aligned.sortedByCoord.dedup.out_plus.bw', OUTDIR=config['OUTDIR'], sample=SAMPLES, REF=["STARsolo_rRNA", "STARsolo"]), # strand-split bigWigs
        # expand('{OUTDIR}/{sample}/{REF}/Aligned.sortedByCoord.dedup.out_merged.bw', OUTDIR=config['OUTDIR'], sample=SAMPLES, REF=["STARsolo_rRNA", "STARsolo"]), #
        expand('{OUTDIR}/{sample}/preTrim_fastqc_{READ}_out', OUTDIR=config['OUTDIR'], sample=SAMPLES, READ=["R1","R2"]), # raw R2 fastQC results
        expand('{OUTDIR}/{sample}/postTrim_fastqc_{READ}_out', OUTDIR=config['OUTDIR'], sample=SAMPLES, READ=["R1","R2"]), # adapter/polyA/ployG-trimmed R1/R2 fastQC results
        expand('{OUTDIR}/{sample}/{REF}/Aligned.sortedByCoord.out.bam.bai', OUTDIR=config['OUTDIR'], sample=SAMPLES, REF=["STARsolo_rRNA", "STARsolo"]), #non-deduplicated .bam; used for saturation estimation
        expand('{OUTDIR}/{sample}/Unmapped_fastqc_out', OUTDIR=config['OUTDIR'], sample=SAMPLES), #fastQC results for unmapped reads
        expand('{OUTDIR}/{sample}/Unmapped.out.mate2_blastResults.txt', OUTDIR=config['OUTDIR'], sample=SAMPLES), # blastn results for unmapped R1 reads non-strand-split bigWigs (for
# fastq preprocessing & QC
include: "rules/1_mergefqs.smk"
include: "rules/1_trimQC.smk"
include: "rules/1_split_bb.smk"

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
# - Dedup .bam file

