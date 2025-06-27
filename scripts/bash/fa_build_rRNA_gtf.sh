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
echo "Pattern for matching: gene_biotype:(${KEYWORDS_PATTERN})"

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
    -v FS=" " \
    -vOFS='' \
    -v pattern="gene_biotype:(${KEYWORDS_PATTERN})" \
    '$5 ~ pattern {
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
        
        split($4, b, ":");
        gene_id = b[2];
        split($6, c, ":");
        transcript_type = c[2];
        split($7, d, ":");
        transcript_name = d[2];
        split($7, e, ":");
        gene_name = e[2];
        split($5, f, ":");
        gene_type = f[2];

        attributes = "gene_id \"" gene_id "\"; transcript_id \"" seqname "\"; gene_type \"" gene_type "\"; gene_name \"" gene_name "\"; transcript_type \"" transcript_type "\"; transcript_name \"" transcript_name "\";";

        print seqname,"\t", source,"\t", feature,"\t", start,"\t", end,"\t", score,"\t", strand,"\t", frame,"\t", attributes;
        count++;
    } END { print "Generated " count " GTF entries" > "/dev/stderr" }' \
| tee >(echo "Compressing and writing GTF output..." > /dev/stderr) \
| gzip \
> ${GTF_rRNA}

# Count GTF entries
echo "Counting generated GTF entries..."
GTF_ENTRIES=$(zcat "${GTF_rRNA}" | wc -l)
echo "Successfully generated ${GTF_ENTRIES} GTF entries"

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