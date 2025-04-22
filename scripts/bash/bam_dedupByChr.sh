#!/bin/bash

# Default values
INBAM=""
CORE=""
OUTBAM=""
TMPDIR=""
CELL_TAG="CB"
UMI_TAG="UB"

# Function to display usage information
_usage() {
cat << EOF
Usage: $0 --bam|-b BAM_FILE --cores|-c CORES --outbam|-o OUTPUT_BAM_FILE --tmpdir|-t TMP_DIRECTORY [--celltag|-ct CELL_TAG] [--umitag|-ut UMI_TAG]

Options:
 --bam|-b BAM_FILE          Path to the input BAM file.
 --cores|-c CORES           Number of cores for parallelization.
 --outbam|-o OUTPUT_BAM_FILE Output BAM file name.
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
        -o|--outbam)
            OUTBAM=$2
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

# Function to log messages
log_message() {
    local message="$1"
    printf "%s - %s\n" "$(date '+%Y-%m-%d %H:%M:%S')" "$message"
}

# Print input parameters
echo "Input parameters:"
echo "  BAM file:          ${INBAM}"
echo "  Cores:             ${CORE}"
echo "  Output BAM file:   ${OUTBAM}"
echo "  Temp. directory:   ${TMPDIR}"
echo "  Cell tag:          ${CELL_TAG}"
echo "  UMI tag:           ${UMI_TAG}"
echo ""

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
    echo "Usage: $0 -b input.bam -c cores -o output.bam -t tmp_directory"
    exit 1
fi

# Check if the output BAM file is provided
if [ -z ${OUTBAM} ]; then
    echo "Error: Output BAM file must be specified using --outbam|-o."
    exit 1
fi

# Ensure the output and temporary directories exist
mkdir -p $(dirname ${OUTBAM})
mkdir -p ${TMPDIR}

# Clean up any leftover files in the temporary directory
log_message "Cleaning up temporary directory: ${TMPDIR}"
rm -rf ${TMPDIR}/*
log_message "Temporary directory cleaned."

PREFIX=$(echo ${INBAM} | rev | cut -d / -f 1 | cut -d . -f 2- | rev)

log_message "Using ${PREFIX} as file prefix..."

# Check for .bam index file
if [ ! -f ${INBAM}.bai ] && [ ! -f ${INBAM}.csi ]; then
    log_message "Indexing input .bam file because you forgot to..."
    samtools index -@ ${CORE} ${INBAM}
fi

# Remove reads that don't have a barcode (CB)
log_message "Removing reads without '${CELL_TAG}' or '${UMI_TAG}' tags..."

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

log_message "Indexing filtered .bam file..."
samtools index -@ ${CORE} ${TMPDIR}/filter.bam

# Check if there are reads remaining in the filtered BAM file
if [ $(samtools view -c ${TMPDIR}/filter.bam) -eq 0 ]; then
    log_message "Error: No reads remaining after filtering. Exiting."
    exit 1
fi

# Extract chromosome names from the filtered BAM file
CHROMOSOMES=$(samtools idxstats ${TMPDIR}/filter.bam | cut -f1)
CHROMOSOMES=$(echo "$CHROMOSOMES" | tr -d '*')

log_message "Found these chromosomes: ${CHROMOSOMES}"

# Loop over each chromosome, split the BAM file, and deduplicate
for CHR in ${CHROMOSOMES}; do
    if [[ "${CHR}" == "*" ]]; then
        log_message "Skipping chromosome '*'."
    else
        log_message "Writing ${TMPDIR}/${CHR}.bam"
        samtools view -b ${TMPDIR}/filter.bam ${CHR} > ${TMPDIR}/${CHR}.bam

        log_message "Deduplicating ${TMPDIR}/${CHR}.bam"
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
        
        log_message "Finished deduplicating ${TMPDIR}/${CHR}.bam"

        rm ${TMPDIR}/${CHR}.bam
    fi
done

# Merge all deduplicated BAM files
log_message "Merging deduplicated BAM files..."
samtools merge -@ ${CORE} ${OUTBAM} ${TMPDIR}/*.dedup.bam

# Index the merged deduplicated BAM file
log_message "Indexing merged deduplicated BAM file..."
samtools index -@ ${CORE} ${OUTBAM}

# Clean up temporary files
if ! rm -r ${TMPDIR}; then
    log_message "Warning: Failed to clean up temporary directory: ${TMPDIR}" >>2
fi

log_message "BAM file deduplicated and split by chromosome into ${OUTBAM}."
