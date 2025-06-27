#!/bin/bash

# Check if the correct number of arguments was provided
if [ "$#" -lt 2 ] || [ "$#" -gt 3 ]; then
    echo "Usage: scripts/bash/extract_rRNA_gtf.sh FASTA_cDNA GTF_rRNA [KEYWORDS]"
    echo "  FASTA_cDNA: Input cDNA FASTA file (compressed or uncompressed)"
    echo "  GTF_rRNA: Output rRNA GTF file (will be gzip compressed)"
    echo "  KEYWORDS: Optional comma-separated list of rRNA keywords (default: rRNA,Mt_rRNA,ribosomal_RNA,5S_rRNA,5.8S_rRNA,18S_rRNA,28S_rRNA,12S_rRNA,16S_rRNA)"
    exit 1
fi

# Assign arguments
FASTA_cDNA=$1
GTF_rRNA=$2

# Set default rRNA keywords or use provided ones
if [ "$#" -eq 3 ]; then
    KEYWORDS=$3
else
    KEYWORDS="rRNA,Mt_rRNA,ribosomal_RNA,5S_rRNA,5.8S_rRNA,18S_rRNA,28S_rRNA,12S_rRNA,16S_rRNA"
fi

# Start timing and logging
START_TIME=$(date +%s)
echo "=== rRNA GTF Generation Started at $(date) ==="
echo "Input file: ${FASTA_cDNA}"
echo "Output file: ${GTF_rRNA}"
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

# Check if headers use pipe-delimited format (like Ensembl) or space-delimited (like GENCODE)
if [[ ${SAMPLE_HEADER} == *"|"* ]]; then
    echo "Detected pipe-delimited format (Ensembl-style)"
    DELIMITER="|"
    FIELD_SEP="|"
    echo "Will search in last field for gene type"
else
    echo "Detected space-delimited format (GENCODE-style)"
    DELIMITER=" "
    FIELD_SEP=" "
    echo "Will search in field 5 for gene_biotype"
fi

# Manually build gtf with this hideous awk code...
# Handle both compressed and uncompressed input files
echo "Building GTF annotations from rRNA sequences..."
if [[ ${FASTA_cDNA} == *.gz ]]; then
    zcat ${FASTA_cDNA}
else
    cat ${FASTA_cDNA}
fi \
| awk \
    -v RS=">" \
    -v FS="${FIELD_SEP}" \
    -vOFS='' \
    -v pattern="(${KEYWORDS_PATTERN})" \
    -v delimiter="${DELIMITER}" \
    'BEGIN { count = 0 }
    NR > 1 {
        # Check if this sequence matches rRNA pattern
        match_found = 0
        if (delimiter == "|") {
            # For pipe-delimited (Ensembl): check last field
            if (NF > 0 && $NF ~ pattern) {
                match_found = 1
            }
        } else {
            # For space-delimited (GENCODE): check field 5 for gene_biotype
            if ($5 ~ ("gene_biotype:" pattern)) {
                match_found = 1
            }
        }
        
        if (match_found) {
            seqname = $1;
            source = "custom";
            feature = "exon";
            start = "1";

            n=split($0, lines, "\n");
            len = 0;
            for (i = 2; i <= n; ++i) {
                len += length(lines[i]);
            }
            end = len;

            score = ".";
            strand = "+";
            frame = ".";
            
            # Handle different header formats for attribute extraction
            if (delimiter == "|") {
                # For pipe-delimited (Ensembl): extract from pipe-separated fields
                # Format: >ENSMUST00000082908.3|ENSMUSG00000064842.3|-|-|Gm26206-201|Gm26206|110|snRNA|
                gene_id = ($2 != "") ? $2 : "unknown"
                transcript_id = ($1 != "") ? $1 : "unknown"
                gene_name = ($6 != "") ? $6 : "unknown"
                transcript_name = ($5 != "") ? $5 : "unknown"
                gene_type = ($NF != "") ? $NF : "unknown"
                transcript_type = gene_type
            } else {
                # For space-delimited (GENCODE): extract from colon-separated values
                split($4, b, ":");
                gene_id = (length(b) > 1) ? b[2] : "unknown";
                split($6, c, ":");
                transcript_type = (length(c) > 1) ? c[2] : "unknown";
                split($7, d, ":");
                transcript_name = (length(d) > 1) ? d[2] : "unknown";
                gene_name = transcript_name;
                split($5, f, ":");
                gene_type = (length(f) > 1) ? f[2] : "unknown";
            }

            attributes = "gene_id \"" gene_id "\"; transcript_id \"" transcript_id "\"; gene_type \"" gene_type "\"; gene_name \"" gene_name "\"; transcript_type \"" transcript_type "\"; transcript_name \"" transcript_name "\";";

            print seqname,"\t", source,"\t", feature,"\t", start,"\t", end,"\t", score,"\t", strand,"\t", frame,"\t", attributes;
            count++;
        }
    } END { 
        print "Generated " count " GTF entries" > "/dev/stderr"
        if (count == 0) {
            print "WARNING: No GTF entries generated!" > "/dev/stderr"
        }
    }' \
| tee >(echo "Compressing and writing GTF output..." > /dev/stderr) \
| gzip \
> ${GTF_rRNA}

# Count GTF entries
echo "Counting generated GTF entries..."
GTF_ENTRIES=$(zcat "${GTF_rRNA}" | wc -l)
echo "Successfully generated ${GTF_ENTRIES} GTF entries"

# Check if output is empty and throw error
if [ "${GTF_ENTRIES}" -eq 0 ]; then
    echo "ERROR: No GTF entries were generated!"
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
    rm -f "${GTF_rRNA}"  # Clean up empty output file
    exit 1
fi

# Calculate runtime and output info
END_TIME=$(date +%s)
RUNTIME=$((END_TIME - START_TIME))
OUTPUT_SIZE=$(du -h "${GTF_rRNA}" | cut -f1)

echo "=== rRNA GTF Generation Summary ==="
echo "Runtime: ${RUNTIME} seconds"
echo "Input sequences: ${TOTAL_SEQS}"
echo "Generated GTF entries: ${GTF_ENTRIES}"
echo "Output file size: ${OUTPUT_SIZE}"
echo "Completed at $(date)"
echo "=== End Summary ==="