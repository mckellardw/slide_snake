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
echo "Pattern for matching: (${KEYWORDS_PATTERN})"

# Analyze header format first
echo "Analyzing header format..."
if [[ ${FASTA_cDNA} == *.gz ]]; then
    SAMPLE_HEADER=$(zcat "${FASTA_cDNA}" | head -1)
else
    SAMPLE_HEADER=$(head -1 "${FASTA_cDNA}")
fi
echo "Sample header: ${SAMPLE_HEADER}"

# Check if headers use pipe-delimited format (GENCODE newer format) or space-delimited (GENCODE older format)
if [[ ${SAMPLE_HEADER} == *"|"* ]]; then
    echo "Detected pipe-delimited format (GENCODE pipe-separated)"
    DELIMITER="|"
    FIELD_SEP="|"
    echo "Will search in last field for gene type"
else
    echo "Detected space-delimited format (GENCODE space-separated)"
    DELIMITER=" "
    FIELD_SEP=" "
    echo "Will search in field 5 for gene_biotype"
fi

# Fields for fasta are different for GENCODE vs Ensembl
# Handle both compressed and uncompressed input files
echo "Extracting rRNA sequences..."
if [[ ${FASTA_cDNA} == *.gz ]]; then
    zcat ${FASTA_cDNA}
else
    cat ${FASTA_cDNA}
fi \
| awk \
    -v RS="\n>" \
    -v FS="${FIELD_SEP}" \
    -v pattern="(${KEYWORDS_PATTERN})" \
    -v delimiter="${DELIMITER}" \
    'BEGIN { count = 0 }
    NR > 1 {
        # For pipe-delimited (GENCODE newer format): check last non-empty field
        if (delimiter == "|") {
            # Remove trailing pipe if present and check the actual last field
            header_line = $0
            gsub(/\|$/, "", header_line)
            split(header_line, fields, "|")
            last_field = fields[length(fields)]
            if (last_field ~ pattern) {
                print ">" $0  # Keep original format
                count++
            }
        }
        # For space-delimited (GENCODE older format): check field 5 for gene_biotype
        else {
            if ($5 ~ ("gene_biotype:" pattern)) {
                print ">" $0
                count++
            }
        }
    }
    END { 
        print "Extracted " count " rRNA sequences" > "/dev/stderr"
        if (count == 0) {
            print "WARNING: No rRNA sequences found!" > "/dev/stderr"
        }
    }' \
| tee >(echo "Compressing and writing output..." > /dev/stderr) \
| gzip \
> ${FASTA_rRNA}

# Count extracted sequences
echo "Counting extracted sequences..."
EXTRACTED_SEQS=$(zcat "${FASTA_rRNA}" | grep -c "^>")
echo "Successfully extracted ${EXTRACTED_SEQS} rRNA sequences"

# Check if output is empty and throw error
if [ "${EXTRACTED_SEQS}" -eq 0 ]; then
    echo "ERROR: No rRNA sequences were extracted!"
    echo "This could be due to:"
    echo "  1. No rRNA sequences in the input file"
    echo "  2. Different header format than expected"
    echo "  3. Keywords don't match the gene types in your file"
    echo ""
    echo "Please check:"
    echo "  - Header format in your input file"
    echo "  - Available gene types in your FASTA headers"
    echo "  - RRNA_KEYWORDS configuration"
    echo ""
    echo "Sample header from your file: ${SAMPLE_HEADER}"
    rm -f "${FASTA_rRNA}"  # Clean up empty output file
    exit 1
fi

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