# fq_readIDs.awk
## Extracts read IDs from a FASTQ - defaults to keeping the 
## first (space-delimited) substring in the header

# Usage:
## zcat file.fq.gz | awk -f fq_readIDs.awk > read_IDs.list

# Extract read ID from header lines (every 4th line starting from line 1)
{if(NR%4==1) split($0,a," "); sub("@", "", a[1]); print a[1]}