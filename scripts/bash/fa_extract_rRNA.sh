#!/bin/bash

# Check if the correct number of arguments was provided
if [ "$#" -lt 2 ] || [ "$#" -gt 3 ]; then
    echo "Usage: scripts/bash/extract_rRNA_fasta.sh FASTA_cDNA FASTA_rRNA [KEYWORDS]"
    echo "  FASTA_cDNA: Input cDNA FASTA file (compressed or uncompressed)"
    echo "  FASTA_rRNA: Output rRNA FASTA file (will be gzip compressed)"
    echo "  KEYWORDS: Optional comma-separated list of rRNA keywords (default: rRNA,Mt_rRNA,ribosomal_RNA,5S_rRNA,5.8S_rRNA,18S_rRNA,28S_rRNA,12S_rRNA,16S_rRNA)"
    exit 1
fi

# Assign arguments
FASTA_cDNA=$1
FASTA_rRNA=$2

# Set default rRNA keywords or use provided ones
if [ "$#" -eq 3 ]; then
    KEYWORDS=$3
else
    KEYWORDS="rRNA,Mt_rRNA,ribosomal_RNA,5S_rRNA,5.8S_rRNA,18S_rRNA,28S_rRNA,12S_rRNA,16S_rRNA"
fi

# Start timing and logging
START_TIME=$(date +%s)
echo "=== rRNA FASTA Extraction Started at $(date) ==="
echo "Input file: ${FASTA_cDNA}"
echo "Output file: ${FASTA_rRNA}"
echo "Keywords: ${KEYWORDS}"

# Check input file exists and get info
if [ ! -f "${FASTA_cDNA}" ]; then
    echo "ERROR: Input file ${FASTA_cDNA} does not exist!"
    exit 1
fi

INPUT_SIZE=$(du -h "${FASTA_cDNA}" | cut -f1)
echo "Input file size: ${INPUT_SIZE}"

# Count total sequences in input
echo "Counting total sequences in input file..."
if [[ ${FASTA_cDNA} == *.gz ]]; then
    TOTAL_SEQS=$(zcat "${FASTA_cDNA}" | grep -c "^>")
else
    TOTAL_SEQS=$(grep -c "^>" "${FASTA_cDNA}")
fi
echo "Total sequences in input: ${TOTAL_SEQS}"

# Create AWK pattern for flexible keyword matching
KEYWORDS_PATTERN=$(echo "${KEYWORDS}" | sed 's/,/|/g')
echo "Pattern for matching: gene_biotype:(${KEYWORDS_PATTERN})"

# Fields for fasta are different for GENCODE
# Handle both compressed and uncompressed input files
echo "Extracting rRNA sequences..."
if [[ ${FASTA_cDNA} == *.gz ]]; then
    zcat ${FASTA_cDNA}
else
    cat ${FASTA_cDNA}
fi \
| awk \
    -v RS="\n>" \
    -v FS=" " \
    -v pattern="gene_biotype:(${KEYWORDS_PATTERN})" \
    '$5 ~ pattern { print ">"$0; count++ } END { print "Extracted " count " rRNA sequences" > "/dev/stderr" }' \
| tee >(echo "Compressing and writing output..." > /dev/stderr) \
| gzip \
> ${FASTA_rRNA}

# Count extracted sequences
echo "Counting extracted sequences..."
EXTRACTED_SEQS=$(zcat "${FASTA_rRNA}" | grep -c "^>")
echo "Successfully extracted ${EXTRACTED_SEQS} rRNA sequences"

# Calculate runtime and output info
END_TIME=$(date +%s)
RUNTIME=$((END_TIME - START_TIME))
OUTPUT_SIZE=$(du -h "${FASTA_rRNA}" | cut -f1)

echo "=== rRNA FASTA Extraction Summary ==="
echo "Runtime: ${RUNTIME} seconds"
echo "Input sequences: ${TOTAL_SEQS}"
echo "Extracted rRNA sequences: ${EXTRACTED_SEQS}"
echo "Extraction rate: $(echo "scale=2; ${EXTRACTED_SEQS}*100/${TOTAL_SEQS}" | bc -l)%"
echo "Output file size: ${OUTPUT_SIZE}"
echo "Completed at $(date)"
echo "=== End Summary ==="