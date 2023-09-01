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

# SAMPLE_SHEET = SAMPLE_SHEET[~SAMPLE_SHEET['sampleID'].str.contains("STO")]
# SAMPLE_SHEET = SAMPLE_SHEET[SAMPLE_SHEET['sampleID'].str.contains("Vis_yPAP_3C")]


SAMPLES = list(SAMPLE_SHEET['sampleID'])

R1_FQS = dict(zip(SAMPLES, list(SAMPLE_SHEET['fastq_R1'])))
R2_FQS = dict(zip(SAMPLES, list(SAMPLE_SHEET['fastq_R2'])))

########################################################################################################
# Executables
########################################################################################################
BWA_EXEC = config['BWA_EXEC']
STAR_EXEC = config['STAR_EXEC']
# KB_EXEC = config['KB_EXEC']
KALLISTO_EXEC = config['KALLISTO_EXEC']
BUST_EXEC = config['BUST_EXEC']
FASTQC_EXEC = config["FASTQC_EXEC"]
CUTADAPT_EXEC = config["CUTADAPT_EXEC"]
SAMTOOLS_EXEC = config["SAMTOOLS_EXEC"]
UMITOOLS_EXEC = config["UMITOOLS_EXEC"]
QUALIMAP_EXEC = config["QUALIMAP_EXEC"]
# MULTIQC_EXEC = config["MULTIQC_EXEC"]
MIRGE_EXEC = config['MIRGE_EXEC']
BOWTIE2_EXEC = config['BOWTIE2_EXEC']
BAM2SPLITBW = config["BAM2SPLITBW"]
FASTX_COLLAPSER = config["FASTX_COLLAPSER"]
BLASTDB = config["BLASTDB"]

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
BB_DICT = {}        # Dictionary of bead barcode maps
SPECIES_DICT = {}   # Dictionary of species listed for mirge3 analysis
for i in range(0,SAMPLE_SHEET.shape[0]):
    tmp_sample = list(SAMPLE_SHEET["sampleID"])[i]
    RECIPE_DICT[tmp_sample] = list(SAMPLE_SHEET["recipe"])[i]
    rRNA_DICT[tmp_sample] = list(SAMPLE_SHEET["STAR_rRNA_ref"])[i]
    REF_DICT[tmp_sample] = list(SAMPLE_SHEET["STAR_ref"])[i]
    GTF_DICT[tmp_sample] = list(SAMPLE_SHEET["genes_gtf"])[i]
    IDX_DICT[tmp_sample] = list(SAMPLE_SHEET["kb_idx"])[i]
    T2G_DICT[tmp_sample] = list(SAMPLE_SHEET["kb_t2g"])[i]
    BB_DICT[tmp_sample] = list(SAMPLE_SHEET["BB_map"])[i]
    SPECIES_DICT[tmp_sample] = list(SAMPLE_SHEET["species"])[i]

########################################################################################################
rule all:
    input:
        # expand( # miRge3.0 pseudobulk analysis
        #     '{OUTDIR}/{sample}/miRge_bulk/annotation.report.html',
        #     OUTDIR=config['OUTDIR'],
        #     sample=SAMPLES
        # ),
        expand( # bowtie2 alignment to small RNA reference(s)
            '{OUTDIR}/{sample}/{SMALL_RNA}/counts.{TYPE}',
            OUTDIR=config['OUTDIR'],
            sample=SAMPLES,
            SMALL_RNA=['piRNA','miRNA'],
            TYPE=["tsv.gz","npz"]
        ),
        expand( # anndata files (with spatial info)
            '{OUTDIR}/{sample}/{ALIGN_OUT}',
            OUTDIR=config['OUTDIR'],
            ALIGN_OUT=['kb/counts_unfiltered/output.h5ad','STARsolo/Solo.out/GeneFull/raw/UniqueAndMultEM.h5ad'],
            sample=SAMPLES
        ),
        expand( #STAR count mats
            '{OUTDIR}/{sample}/{REF}/Solo.out/GeneFull/raw/matrix.mtx.gz',
            OUTDIR=config['OUTDIR'],
            sample=SAMPLES,
            REF=["STARsolo_rRNA", "STARsolo"]
        ),
        expand( # kallisto/bustools count mats
            '{OUTDIR}/{sample}/kb/counts_unfiltered/output.mtx.gz',
            OUTDIR=config['OUTDIR'],
            sample=SAMPLES
        ),
        # expand( #TODO- REF=["STARsolo_rRNA", "STARsolo"]), # umi_tools deduplicated .bam
        #     '{OUTDIR}/{sample}/STARsolo/Aligned.sortedByCoord.dedup.out.bam.bai',
        #     OUTDIR=config['OUTDIR'],
        #     sample=SAMPLES
        # ),
        expand(  # alignment QC with qualimap | requires deduped input!
            '{OUTDIR}/{sample}/qualimap/{FILE}',
            OUTDIR=config['OUTDIR'],
            sample=SAMPLES,
            FILE=["qualimapReport.html","rnaseq_qc_result.csv"]
        ),
        # expand( # strand-split, umi_tools deduplicated .bam
        #     '{OUTDIR}/{sample}/STARsolo/Aligned.sortedByCoord.dedup.out.{STRAND}.bam.bai',
        #     OUTDIR=config['OUTDIR'],
        #     sample=SAMPLES,
        #     STRAND=["fwd", "rev"]
        # ),
        # expand( #non-deduplicated .bam
        #     '{OUTDIR}/{sample}/{REF}/Aligned.sortedByCoord.out.bam.bai',
        #     OUTDIR=config['OUTDIR'],
        #     sample=SAMPLES,
        #     REF=["STARsolo_rRNA", "STARsolo"]
        # ),
        expand( #fastQC results for unmapped reads
            '{OUTDIR}/{sample}/fastqc_unmapped',
            OUTDIR=config['OUTDIR'],
            sample=SAMPLES
        ),
        # expand( # blastn results for unmapped R2 reads
        #     '{OUTDIR}/{sample}/Unmapped.out.mate2_blastResults.txt',
        #     OUTDIR=config['OUTDIR'],
        #     sample=SAMPLES
        # ),
        expand( # fastQC results
            '{OUTDIR}/{sample}/fastqc_{TRIM}_{READ}',
            OUTDIR=config['OUTDIR'],
            sample=SAMPLES,
            READ=["R1","R2"],
            TRIM = ["preTrim","postTrim"]
        )
        

# fastq preprocessing & QC
include: "rules/1a_mergefqs.smk"
include: "rules/1b_trimQC.smk"
include: "rules/1c_split_bb.smk"

# STAR alignment, QC, and post-processing
include: "rules/2a_rRNA_filter.smk"
include: "rules/2b_star_align.smk"
include: "rules/2c_star_unmapped.smk"
include: "rules/2d_star_qualimap.smk"
include: "rules/2e_star_dedup.smk"

# kallisto/bustools alignment
include: "rules/3a_kallisto_align.smk"
include: "rules/3b_kallisto_pseudobam.smk"

# small RNA stuff
include: "rules/4_mirge.smk"
include: "rules/4_piRNA_bowtie2.smk"
include: "rules/4_miRNA_bowtie2.smk"

# scanpy stuff
include: "rules/5a_scanpy_init.smk"

# Post-processing, prep for downstream analyses
#TODO:
# - Add proper saturation estimation (not just saturation based on reads aligning to known genes)
