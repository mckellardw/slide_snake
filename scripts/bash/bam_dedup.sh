#!/usr/bin/bash
# bam_dedup.sh - a bash script that deduplicates .bam files, aligned with STARsolo
#                - deduplicates chr-by-chr, to reduce run time and memory requirements
# Usage:
# bash bam_dedup.sh --bam|-b /path/to/input.bam --whitelist|-w /path/to/whitelist --cores|-c num_cores --outbam|-o /path/to/output.bam --tmpdir|-t /path/to/tmp/directory [--celltag|-ct CELL_TAG] [--umitag|-ut UMI_TAG]

# Default values
INBAM=""
WHITELIST=""
CORE=""
OUTBAM=""
TMPDIR=""
CELL_TAG="CB"
UMI_TAG="UB"

# Function to display usage information
_usage() {
cat << EOF
Usage: $0 --bam|-b BAM_FILE --whitelist|-w WHITELIST_FILE --cores|-c CORES --outbam|-o OUTPUT_BAM --tmpdir|-t TMP_DIRECTORY [--celltag|-ct CELL_TAG] [--umitag|-ut UMI_TAG]

Options:
 --bam|-b BAM_FILE          Path to the input BAM file.
 --whitelist|-w WHITELIST_FILE Path to the barcode whitelist.
 --cores|-c CORES           Number of cores for parallelization.
 --outbam|-o OUTPUT_BAM     Path to the output deduplicated BAM file.
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
        -w|--whitelist)
            WHITELIST=$2
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

# Print input parameters
echo "Input parameters:"
echo "  BAM file:          ${INBAM}"
echo "  Whitelist file:    ${WHITELIST}"
echo "  Cores:             ${CORE}"
echo "  Output BAM:        ${OUTBAM}"
echo "  Temp. directory:   ${TMPDIR}"
echo "  Cell tag:          ${CELL_TAG}"
echo "  UMI tag:           ${UMI_TAG}"
echo ""

# Check if required parameters are provided
if [ -z ${INBAM} ] || [ -z ${WHITELIST} ] || [ -z ${CORE} ] || [ -z ${OUTBAM} ] || [ -z ${TMPDIR} ]; then
    echo "Error: Missing required arguments."
    _usage
    exit 1
fi

# Ensure the temporary directory exists
mkdir -p ${TMPDIR}

PREFIX=$(echo ${INBAM} | rev | cut -d / -f 1 | cut -d . -f 2- | rev)

echo "Using ${PREFIX} as file prefix..."

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
    ${INBAM} \
| grep -v "${CELL_TAG}:Z:-" \
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
    --temp-dir=${TMPDIR} \
    --multimapping-detection-method=NH \
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
