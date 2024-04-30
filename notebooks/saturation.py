
import pandas as pd
import os
import yaml


print("Loading Parameters")

with open("../../config/config.yaml") as stream:
    try:
        config=yaml.safe_load(stream)
    except yaml.YAMLError as exc:
        print(exc)

print(".")
#Parameters
OUTDIR='/gpfs/commons/home/ecordina/uST/github/slide_snake/out'
SAMPLE="decoderseq"
RECIPE='decoderseq'
RECIPE_DICT=pd.read_csv("../../resources/eman_recipe_sheet.csv",index_col=0)

#input/Temp output
print("..")
BB_3=OUTDIR+"/"+SAMPLE+"/bb/whitelist_underscore.txt"
BAM = OUTDIR+"/"+SAMPLE+"/STARsolo/short_read/"+RECIPE+"/Aligned.sortedByCoord.out.bam"
DOWNSAMPLE1 = OUTDIR+"/"+SAMPLE+"/STARsolo/short_read/"+RECIPE+"/Downsample1Aligned.sortedByCoord.out.bam"
DOWNSAMPLE2 = OUTDIR+"/"+SAMPLE+"/STARsolo/short_read/"+RECIPE+"/Downsample2Aligned.sortedByCoord.out.bam"
BAMDEDUP=OUTDIR+"/"+SAMPLE+"/STARsolo/short_read/"+RECIPE+"/Aligned.sortedByCoord.dedup.out.bam"
DOWNSAMPLE1DEDUP = OUTDIR+"/"+SAMPLE+"/STARsolo/short_read/"+RECIPE+"/Downsample1Aligned.sortedByCoord.dedup.out.bam"
DOWNSAMPLE2DEDUP = OUTDIR+"/"+SAMPLE+"/STARsolo/short_read/"+RECIPE+"/Downsample2Aligned.sortedByCoord.dedup.out.bam"
#Temps
TMP_DIR = OUTDIR+"/"+SAMPLE+"/tmp"
log=OUTDIR+"/"+SAMPLE+"/STARsolo/short_read/"+RECIPE+"/dedup.log"
#Output
SAT = OUTDIR+"/"+SAMPLE+"/STARsolo/short_read/"+RECIPE+"/saturation.csv"

print("...")
threads=config["CORES"]
tmp_recipe = RECIPE_DICT.loc[SAMPLE]
INPUTDIR="/".join(BAM.split("/")[:-1])+"/"
subsample=[42.5,42.3]

print("Parameters loaded!")
print("Starting SubSampling")
if not os.path.isfile(DOWNSAMPLE1):
    os.system(
        f"""
        samtools view \
            -s {subsample[0]} \
            --threads {threads} \
            --bam {BAM} \
            --output  {DOWNSAMPLE1}
            """
        )
    print("First Subsampling done!")
if not os.path.isfile(DOWNSAMPLE2):
    os.system(
        f"""
        samtools view \
            -s {subsample[1]} \
            --threads {threads} \
            --bam {BAM} \
            --output  {DOWNSAMPLE2}
            """
        )
print("Subsampling Finished!")
print("Starting deduplication")

if not os.path.isfile(DOWNSAMPLE1DEDUP):
    os.system(
        f"""
        bash ../bash/split_dedup.sh \
            {DOWNSAMPLE1} \
            {BB_3} \
            {threads} \
            {DOWNSAMPLE1DEDUP} \
            $(dirname {DOWNSAMPLE1DEDUP})/tmp/dedup \
            | tee {log}
        """
    )
print("first Deduplication Done!")
if not os.path.isfile(DOWNSAMPLE2DEDUP):
    os.system(
        f"""
        bash ../bash/split_dedup.sh \
            {DOWNSAMPLE2} \
            {BB_3} \
            {threads} \
            {DOWNSAMPLE2DEDUP} \
            $(dirname {DOWNSAMPLE2DEDUP})/tmp/dedup \
            | tee {log}
        """
        )
print("Deduplication Done!")
print("Computing number of reads")
if not os.path.isfile("tot_nreads.txt"):
    os.system(
        f"""
        samtools view -c {BAM} > tot_nreads.txt
        """
        )
if not os.path.isfile("nreads.txt"):
    os.system(
        f"""
        samtools view -c {BAMDEDUP} > nreads.txt
        """
        )
if not os.path.isfile("nreads_ds1.txt"):
    os.system(
        f"""
        samtools view -c {DOWNSAMPLE1DEDUP} > nreads_ds1.txt
        """
        )
if not os.path.isfile("nreads_ds2.txt"):
    os.system(
        f"""
        samtools view -c {DOWNSAMPLE2DEDUP} > nreads_ds2.txt
        """
        )
print("Computing Saturation!")
with open("tot_nreads.txt", "r") as f:
    tot_reads=int(f.read())
dedup_reads=[]
for file in ["nreads.txt","nreads_ds1.txt","nreads_ds2.txt"]:
    with open(file, "r") as f:
        dedup_reads.append(1-int(f.read())/tot_reads)
downsampling=[1,0.5,0.3]
dataFrame=pd.DataFrame(columns=['No Downsample',"Downsample1","Downsample2"],index=["Downsample Value","Saturation Value"],data=[downsampling,dedup_reads])
dataFrame.to_csv(SAT)
print("Done!")
print(f"Saturation levels saved at {SAT}")
