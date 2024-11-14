#!/usr/bin/bash
# bam_dedup.sh - a bash script that deduplicates .bam files, aligned with STARsolo
#                - deduplicates chr-by-chr, to reduce run time and memory requirements
# Usage:
# bash bam_dedup.sh /path/to/input.bam /path/to/whitelist num_cores /path/to/output/directory /path/to/tmp/directory [cell_tag] [umi_tag]

INBAM=$1 # path to .bam file (sorted & indexed already!)
WHITELIST=$2 # path to barcode whitelist
CORE=$3 # number of cores for parallelization
OUTBAM=$4 # output/deduped bam path
TMPDIR=$5
CELL_TAG=${6:-CB} # default to CB if not provided
UMI_TAG=${7:-UB}  # default to UB if not provided

mkdir -p ${TMPDIR}
cd ${TMPDIR}

PREFIX=$(echo ${INBAM} | rev | cut -d / -f 1 | cut -d . -f 2- | rev)

echo "Using ${PREFIX} as file prefix..."

# Initialize log output
echo ".bam file location: ${INBAM}"
echo "Max cores:          ${CORE}"
echo "Output location:    ${OUTBAM}"
echo

# Check params...
if [ ! -f ${INBAM} ]; then
    echo "Can't find ${INBAM}"
fi

# Check for .bam index file
if [ ! -f ${INBAM}.bai ]; then
    echo "Indexing input .bam file because you forgot to..."
    samtools index -@ ${CORE} ${INBAM}
fi

# Remove reads that don't have a barcode (CB)
echo "Removing reads without '${CELL_TAG}' or '${UMI_TAG}' tags..."
date

samtools view \
    -h \
    -@ ${CORE} \
    -q 1 \
    --tag-file ${CELL_TAG}:${WHITELIST} \
    ${INBAM} \
| grep -v "${UMI_TAG}:Z:-" \
| samtools view -bS \
> ${TMPDIR}/filter.bam

echo

# Index split .bam's
echo "Indexing filtered .bam file..."
date

samtools index -@ ${CORE} ${TMPDIR}/filter.bam

echo

# Dedup resulting bams one-by-one
echo "Deduplicating filtered .bam file..."
date

umi_tools dedup \
    -I ${TMPDIR}/filter.bam \
    --extract-umi-method=tag \
    --umi-tag=${UMI_TAG} \
    --cell-tag=${CELL_TAG} \
    --method=unique \
    --per-cell \
    --unmapped-reads=discard \
    -S ${TMPDIR}/dedup.bam

echo

# Merge, sort, and index dedup'ed .bams
echo "Sorting and indexing deduplicated .bam file..."
date

samtools sort ${TMPDIR}/dedup.bam > ${OUTBAM}
samtools index -@ ${CORE} ${OUTBAM} && rm -r ${TMPDIR}
