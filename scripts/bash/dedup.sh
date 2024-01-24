#!/usr/bin/bash
# split_dedup.sh - a bash script that deduplicates .bam files, aligned with STARsolo
#                - deduplicates chr-by-chr, to reduce run time and memory requirements
# Usage:
# bash chr_split_dedup.sh /path/to/input.bam /path/to/whitelist num_cores /path/to/output/directory /path/to/tmp/directory file.LOG

# bash /home/dwm269/DWM_utils/seq_utils/chr_split_dedup.sh ds.bam /workdir/dwm269/totalRNA/STRS-HD/data/STARsolo/SH2A/tmp/whitelist.txt 24 ./OUTPUT_DEDUP DEDUP.LOG

# chr_split_dedup (v1.0) - input a .bam file, output the deduplicated .bam file in OUTDIR

#TODO: python script instead of bash to make parallelization a bit easier?

INBAM=$1 # path to .bam file (sorted & indexed already!)
BB=$2 # path to barcode whitelist
CORE=$3 # number of cores for parallelization
OUTBAM=$4 # output/deduped bam path
TMPDIR=$5
LOG=$6


mkdir -p ${TMPDIR}
cd ${TMPDIR}

echo "Log info can be found in "${LOG}

PREFIX=`echo ${INBAM} | rev | cut -d / -f 1 | cut -d . -f 2- | rev`

echo "Using "${PREFIX}" as file prefix..." 

# Initialize log file
echo ".bam file location:    "${INBAM}  > ${LOG}
echo "Max cores:             "${CORE}   >> ${LOG}
echo "Output location:       "${OUTBAM} >> ${LOG}
echo >> ${LOG}
echo >> ${LOG}

# Check params...
if [ ! -f ${INBAM} ]
then
    echo "Can't find ${INBAM}"
fi

# Check for .bam index file
if [ ! -f ${INBAM}.bai ]
then
    echo "Indexing input .bam file because you forgot to..." >> ${LOG}
    samtools index -@ ${CORE} ${INBAM}
fi

# Remove reads that don't have a barcode (CB)
echo "Removing reads without 'CB' or 'UB' tags..." >> ${LOG}
date >> ${LOG}

samtools view \
-h \
-@ ${CORE} \
-q 1 \
--tag-file CB:${BB} \
${INBAM} \
| grep -v "UB:Z:-" \
| samtools view -bS \
> ${TMPDIR}/filter.bam \
| tee -a {LOG}

echo >> ${LOG}


# index split .bam's
echo "Indexing filtered .bam file..." >> ${LOG}
date >> ${LOG}

samtools index \
-@ ${CORE} \
${TMPDIR}/filter.bam

echo >> ${LOG}

# Dedup resulting bams one-by-one
echo "Deduplicating filtered .bam file..." >> ${LOG}
date >> ${LOG}


umi_tools dedup \
-I ${TMPDIR}/filter.bam \
--extract-umi-method=tag \
--umi-tag=UB \
--cell-tag=CB \
--method=unique \
--per-cell \
--unmapped-reads=discard \
-S ${TMPDIR}/dedup.bam \
| tee -a {LOG}
# --output-stats={params.OUTPUT_PREFIX} \

echo >> {LOG}

# merge, sort, and index dedup'ed .bams
echo "Sorting and indexing deduplicated .bam file..." >> {LOG}
date >> {LOG}

samtools sort ${TMPDIR}/dedup.bam > ${OUTBAM}
samtools index -@ ${CORE} ${OUTBAM} && rm -r ${TMPDIR}
