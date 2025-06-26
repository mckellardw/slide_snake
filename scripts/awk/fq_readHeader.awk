# fq_readHeader.awk
## Extract complete read headers from a FASTQ file (without @ symbol)

# Usage:
## zcat file.fq.gz | awk -f fq_readHeader.awk > headers.txt

# Extract header lines (every 4th line starting from line 1) and remove @ symbol
{if(NR%4==1) print substr($0,2)}