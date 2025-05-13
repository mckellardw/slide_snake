#!/bin/bash

#SBATCH --job-name=pathodbit
#SBATCH --output=.slurm/%j.out
#SBATCH --partition=pe2
#SBATCH --cpus-per-task=24
#SBATCH --mem=250G
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=dmckellar@nygenome.org
#SBATCH --time=02-00:00:00

# Preprint - 
# Paper - https://www.cell.com/cell/fulltext/S0092-8674(24)01019-5
# Code - https://github.com/Zhiliang-Bai/Patho-DBiT
# GEO accession - GSE274641

## conda env:
# mamba activate fastq-dl

DATADIR="/path/to/data/patho_dbit"
mkdir -p ${DATADIR}
cd ${DATADIR}

ACC_LIST="/path/to/slide_snake/data/test_pathodbit/SRR_Acc_List.txt"

NTHREADS=24

## load in SRR IDs
readarray -t SRR < ${ACC_LIST}

for i in "${SRR[@]}"
do
  echo ${i}
  fastq-dl \
    --accession ${i} \
    --outdir ${DATADIR} \
    --prefix ${i} \
    --cpus ${NTHREADS} \
    --group-by-sample
done


## w/ prefetch & parallel-fastq-dump
# for i in "${SRR[@]}"
# do
#   echo ${i}
#   prefetch \
#     --verify yes \
#     --max-size 9999999999 \
#     --output-directory ${DATADIR} \
#     ${i}

#   parallel-fastq-dump \
#     --sra-id ${i}.sra \
#     --threads ${NTHREADS} \
#     --outdir ${DATADIR} \
#     --tmpdir ${DATADIR} \
#     --split-files \
#     --gzip && \
#   rm ${i}.sra
# done
