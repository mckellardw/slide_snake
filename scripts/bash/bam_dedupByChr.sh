#!/bin/bash

# Default values
INBAM=""
CORE=""
OUTDIR=""
TMPDIR=""
CELL_TAG="CB"
UMI_TAG="UB"

# Function to display usage information
_usage() {
cat << EOF
Usage: $0 --bam|-b BAM_FILE --cores|-c CORES --outdir|-o OUTPUT_DIRECTORY --tmpdir|-t TMP_DIRECTORY [--celltag|-ct CELL_TAG] [--umitag|-ut UMI_TAG]

Options:
 --bam|-b BAM_FILE          Path to the input BAM file.
 --cores|-c CORES           Number of cores for parallelization.
 --outdir|-o OUTPUT_DIRECTORY Output directory for split BAM files.
 --tmpdir|-t TMP_DIRECTORY  Temporary directory for intermediate files.
 --celltag|-ct CELL_TAG     Cell tag (default: CB).
 --umitag|-ut UMI_TAG       UMI tag (default: UB).
EOF
}

# Parse command-line arguments
while [[ "$#" -gt 0 ]]; do
    case $1 in
        -b|--bam)
            INBAM=$2
            shift 2
            ;;
        -c|--cores)
            CORE=$2
            shift 2
            ;;
        -o|--outdir)
            OUTDIR=$2
            shift 2
            ;;
        -t|--tmpdir)
            TMPDIR=$2
            shift 2
            ;;
        -ct|--celltag)
            CELL_TAG=$2
            shift 2
            ;;
        -ut|--umitag)
            UMI_TAG=$2
            shift 2
            ;;
        *)
            echo "Unknown argument: $1" >&2
            _usage >&2
            exit 1
            ;;
    esac
done

# Check if samtools and umi_tools are installed
if ! command -v samtools &> /dev/null; then
    echo "samtools could not be found, please install it first"
    exit 1
fi

if ! command -v umi_tools &> /dev/null; then
    echo "umi_tools could not be found, please install it first"
    exit 1
fi

# Check if the BAM file is provided
if [ -z ${INBAM} ]; then
    echo "Usage: $0 -b input.bam -c cores -o output_directory -t tmp_directory"
    exit 1
fi

# Set default output directory to the same directory as the BAM file if not specified
if [ -z ${OUTDIR} ]; then
    OUTDIR=$(dirname ${INBAM})
fi

# Ensure the output and temporary directories exist
mkdir -p ${OUTDIR}
mkdir -p ${TMPDIR}

PREFIX=$(echo ${INBAM} | rev | cut -d / -f 1 | cut -d . -f 2- | rev)

echo "Using ${PREFIX} as file prefix..."

# Check for .bam index file
if [ ! -f ${INBAM}.bai ] && [ ! -f ${INBAM}.csi ]; then
    echo "Indexing input .bam file because you forgot to..."
    samtools index -@ ${CORE} ${INBAM}
fi

# Remove reads that don't have a barcode (CB)
echo "Removing reads without '${CELL_TAG}' or '${UMI_TAG}' tags..."
date

# Remove unmapped (q<1) and untagged (CELL_TAG or UMI_TAG) reads
samtools view \
    -h \
    -@ ${CORE} \
    -q 1 \
    ${INBAM} \
| grep -v "${CELL_TAG}:Z:-" \
| grep -v "${UMI_TAG}:Z:-" \
| samtools view -bS \
> ${TMPDIR}/filter.bam

echo

# Index filtered .bam
echo "Indexing filtered .bam file..."
date

samtools index -@ ${CORE} ${TMPDIR}/filter.bam

# Check if there are reads remaining in the filtered BAM file
if [ $(samtools view -c ${TMPDIR}/filter.bam) -eq 0 ]; then
    echo "Error: No reads remaining after filtering. Exiting."
    exit 1
fi

echo ""

# Extract chromosome names from the filtered BAM file
CHROMOSOMES=$(samtools idxstats ${TMPDIR}/filter.bam | cut -f1)

# Remove asterisk character from CHROMOSOMES
CHROMOSOMES=$(echo "$CHROMOSOMES" | tr -d '*')

echo "Found these chromosomes: "
echo "${CHROMOSOMES}"
echo ""

# Loop over each chromosome, split the BAM file, and deduplicate
for CHR in ${CHROMOSOMES}; do
    if [[ "${CHR}" == "*" ]]; then
        echo "Skipping chromosome '*'."
    else
        echo "    Writing ${TMPDIR}/${CHR}.bam"
        samtools view -b ${TMPDIR}/filter.bam ${CHR} > ${TMPDIR}/${CHR}.bam

        echo "    Deduplicating ${TMPDIR}/${CHR}.bam"
        samtools index -@ ${CORE} ${TMPDIR}/${CHR}.bam  
        umi_tools dedup \
            -I ${TMPDIR}/${CHR}.bam \
            --temp-dir=${TMPDIR} \
            --multimapping-detection-method=NH \
            --extract-umi-method=tag \
            --umi-tag=${UMI_TAG} \
            --cell-tag=${CELL_TAG} \
            --method=unique \
            --per-cell \
            --unmapped-reads=discard \
            -S ${TMPDIR}/${CHR}.dedup.bam
        
        echo ""

        rm ${TMPDIR}/${CHR}.bam
    fi  # Fix syntax error by adding this line
done

# Merge all deduplicated BAM files
echo "Merging deduplicated BAM files..."
samtools merge -@ ${CORE} ${OUTDIR}/${PREFIX}.dedup.bam ${TMPDIR}/*.dedup.bam

# Index the merged deduplicated BAM file
echo "Indexing merged deduplicated BAM file..."
samtools index -@ ${CORE} ${OUTDIR}/${PREFIX}.dedup.bam

# Clean up temporary files
# rm -r ${TMPDIR}

echo "BAM file deduplicated and split by chromosome into ${OUTDIR}."
