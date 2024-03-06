#!/bin/bash

# Default values
BAM_FILE=""
OUTDIR=""

# Function to display usage information
_usage() {
cat << EOF
Usage: $0 --bam|-b BAM_FILE --outdir|-o OUTPUT_DIRECTORY

Options:
 --bam|-b BAM_FILE    Path to the input BAM file.
 --outdir|-o OUTPUT_DIRECTORY Output directory for split BAM files. Defaults to the same directory as the input BAM file.
EOF
}

# Parse command-line arguments
while [[ "$#" -gt 0 ]]; do
    case $1 in
        -b|--bam)
            BAM_FILE=$2
            shift 2
            ;;
        -o|--outdir)
            OUTDIR=$2
            shift 2
            ;;
        *)
            echo "Unknown argument: $1" >&2
            _usage >&2
            exit 1
            ;;
    esac
done

# Check if samtools is installed
if ! command -v samtools &> /dev/null; then
    echo "samtools could not be found, please install it first"
    exit 1
fi

# Check if the BAM file is provided
if [ -z ${BAM_FILE} ]; then
    echo "Usage: $0 -b input.bam [-o output_directory]"
    exit 1
fi

# Set default output directory to the same directory as the BAM file if not specified
if [ -z ${OUTDIR} ]; then
    OUTDIR=$(dirname ${BAM_FILE})
fi

# Ensure the output directory exists
mkdir -p ${OUTDIR}

# Check if the BAM file is indexed
if [ ! -f ${BAM_FILE}.bai ] && [ ! -f ${BAM_FILE}.csi ]; then
    echo "Indexing BAM file..."
    samtools index "$BAM_FILE"
fi

# Extract chromosome names from the BAM file
CHROMOSOMES=$(samtools idxstats ${BAM_FILE} | cut -f1)

# Remove asterisk character from CHROMOSOMES
CHROMOSOMES=$(echo "$CHROMOSOMES" | tr -d '*')

echo "Found these chromosomes: "
echo "${CHROMOSOMES}"
echo ""

# Loop over each chromosome and split the BAM file
for CHR in ${CHROMOSOMES}; do
    if [[ "${CHR}" == "*" ]]; then
        echo "Skipping chromosome '*'."
    else
        # echo "$CHR"
        echo "    Writing ${OUTDIR}/${CHR}.bam"
        samtools view -b ${BAM_FILE} ${CHR} > ${OUTDIR}/${CHR}.bam
    fi
done

# Write unmapped reads to "unmapped.bam"
echo "Writing unmapped reads to 'unmapped.bam'..."
samtools view -b -f 4 ${BAM_FILE} > ${OUTDIR}/unmapped.bam

echo "BAM file split by chromosome into ${OUTDIR}."
