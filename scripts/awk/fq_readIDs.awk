# fq_readIDs.awk
##  Extracts read IDs from a fastq - defaults to keeping the 
##     first (space-delimited) substring in the header

# Usage:
## zcat file.fq.gz | awk  -f fq_readIDs.awk > read_IDs.list

{if(NR%4==1) split($0,a," "); sub("@", "", a[1]); print a[1]}