#!/bin/bash

START_TIME=$(date +%s)

# Usage: ./fa_extract_rRNA.sh <input_cDNA_fasta> <output_rRNA_fasta.gz> [KEYWORDS]
# Example: ./fa_extract_rRNA.sh gencode.vM31.transcripts.fa.gz rRNA.fa.gz

if [ "$#" -lt 2 ] || [ "$#" -gt 3 ]; then
    echo "Usage: $0 <input_cDNA_fasta> <output_rRNA_fasta.gz> [KEYWORDS]"
    echo "  <input_cDNA_fasta>: GENCODE cDNA FASTA file (.fa or .fa.gz)"
    echo "  <output_rRNA_fasta.gz>: Output gzipped FASTA file"
    echo "  [KEYWORDS]: Optional, comma-separated rRNA keywords (default: rRNA,Mt_rRNA,ribosomal_RNA,5S_rRNA,5.8S_rRNA,18S_rRNA,28S_rRNA,12S_rRNA,16S_rRNA)"
    exit 1
fi

IN_FASTA="$1"
OUT_FASTA="$2"
KEYWORDS="${3:-rRNA,Mt_rRNA,ribosomal_RNA,5S_rRNA,5.8S_rRNA,18S_rRNA,28S_rRNA,12S_rRNA,16S_rRNA}"

# Log input information
echo "Input file:   $IN_FASTA"
echo "Output file:  $OUT_FASTA"
echo "Keywords:     $KEYWORDS"
echo ""

# Check if input file exists
if [[ ! -f "$IN_FASTA" ]]; then
    echo "Error: Input file '$IN_FASTA' does not exist."
    exit 1
fi

# Check if output directory exists, create if not
OUT_DIR=$(dirname "$OUT_FASTA")
if [[ ! -d "$OUT_DIR" ]]; then
    echo "Output directory '$OUT_DIR' does not exist. Creating it."
    mkdir -p "$OUT_DIR"
    if [[ $? -ne 0 ]]; then
        echo "Error: Failed to create output directory '$OUT_DIR'."
        exit 1
    fi
fi

# Check if input file is gzipped
if [[ "$IN_FASTA" == *.gz ]]; then
    SAMPLE_HEADER=$(zcat "$IN_FASTA" | head -1)
else
    SAMPLE_HEADER=$(head -1 "$IN_FASTA")
fi
echo "Sample header: $SAMPLE_HEADER"

# Prepare input stream
if [[ "$IN_FASTA" == *.gz ]]; then
    IN_CMD="zcat \"$IN_FASTA\""
else
    IN_CMD="cat \"$IN_FASTA\""
fi

awk_script='
BEGIN {
    n = split(keywords, kw, ",");
    keep = 0;
}
{
    if ($0 ~ /^>/) {
        keep = 0;
        for (i = 1; i <= n; i++) {
            if (index($0, kw[i])) {
                keep = 1;
                break;
            }
        }
    }
    if (keep) print $0;
}
'

eval $IN_CMD | awk -v keywords="$KEYWORDS" "$awk_script" | gzip > "$OUT_FASTA"

COUNT=$(zcat "$OUT_FASTA" | grep -c "^>")
echo "Extracted $COUNT rRNA sequences to $OUT_FASTA"

END_TIME=$(date +%s)
RUNTIME=$((END_TIME - START_TIME))
echo "Runtime: ${RUNTIME} seconds"