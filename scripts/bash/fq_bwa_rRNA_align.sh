#!/bin/bash

# Parse input flags
while getopts "r:a:s:n:f:q:t:l:" opt; do
    case "$opt" in
        r) R2_FQ="$OPTARG";;
        a) BAM1="$OPTARG";;
        s) BAM2="$OPTARG";;
        n) NO_RIBO_FQ="$OPTARG";;
        f) BWA_REF="$OPTARG";;
        q) MIN_MAPQ="$OPTARG";;
        t) THREADS="$OPTARG";;
        *) echo "Usage: $0 -r <R2_FQ> -a <BAM1> -s <BAM2> -n <NO_RIBO_FQ> -f <BWA_REF> -q <MIN_MAPQ> -t <THREADS>"
           exit 1;;
    esac
done

MIN_NM=5

# Log passed parameters for debugging
echo "fq_bwa_rRNA_align:"
echo "Input fastq:             ${R2_FQ}"
echo "Aligned .bam:            ${BAM1}"
echo "Aligned & Sorte .bam:    ${BAM2}"
echo "Output .fastq (no rRNA): ${NO_RIBO_FQ}"
echo "rRNA reference:          ${BWA_REF}"
echo "Minimum MAPQ to remove:  ${MIN_MAPQ}"
echo "Minimum NM to remove:    ${MIN_NM}"
echo "Num threads:             ${THREADS}"
echo " "

# Check for missing arguments
if [ -z "${R2_FQ}" ] || [ -z "${BAM1}" ] || [ -z "${BAM2}" ] || [ -z "${NO_RIBO_FQ}" ] || [ -z "${BWA_REF}" ] || [ -z "${MIN_MAPQ}" ] || [ -z "${THREADS}" ]; then
    echo "Usage: $0 -r <R2_FQ> -a <BAM1> -s <BAM2> -n <NO_RIBO_FQ> -f <BWA_REF> -q <MIN_MAPQ> -t <THREADS>"
    exit 1
fi

# Create output directory if not exists
mkdir -p "$(dirname "${BAM1}")"

# Run bwa-mem2
echo $(date "+%Y-%m-%d %H:%M:%S")": aligning to rRNA reference..."
bwa-mem2 mem -t "${THREADS}" "${BWA_REF}" "${R2_FQ}" \
1> "${BAM1}"

# Sort BAM file
echo $(date "+%Y-%m-%d %H:%M:%S")": sorting output..."
samtools sort -@ "${THREADS}" -O BAM -o "${BAM2}" "${BAM1}"

# Index the sorted BAM file
echo $(date "+%Y-%m-%d %H:%M:%S")": Indexing bam..."
samtools index "${BAM2}"

# Tag	Description
# AS:i:	Alignment score (higher is better)
# XS:i:	Suboptimal alignment score (for secondary alignments)
# NM:i:	Edit distance (number of mismatches and indels)
# MD:Z:	Mismatch string (indicates mismatches against the reference)

# Check if there are aligned reads in BAM2
ALIGNED_COUNT=$(samtools view -c -F 4 "${BAM2}")
if [ "${ALIGNED_COUNT}" -eq 0 ]; then
    echo $(date "+%Y-%m-%d %H:%M:%S")": No aligned reads found. Creating empty .fastq file."
    > "${NO_RIBO_FQ}"
else
    echo $(date "+%Y-%m-%d %H:%M:%S")": Filtering aligned reads..."
    # Filter high edit-distance ("NM">) and low MAPQ reads and convert BAM to fastq
    samtools view -h "${BAM2}" | \
        awk -v tag=NM -v threshold="${MIN_NM}" -f scripts/awk/bam_highPassTagFilter.awk | \
        awk -v quality="${MIN_MAPQ}" -f scripts/awk/bam_lowPassMAPQFilter.awk | \
        samtools fastq - > "${NO_RIBO_FQ}"
fi

echo $(date "+%Y-%m-%d %H:%M:%S")": Done!"