# bam_shortPassReadFilter.awk
## Filter reads in a BAM file by length (short pass filter - keeps reads shorter than threshold)

# Usage:
## samtools view -h input.bam | awk -v max_length=100 -f bam_shortPassReadFilter.awk > output.sam

BEGIN { 
    FS = OFS = "\t" 
}

# Keep header lines and reads shorter than max_length
(length($10) < max_length || $1 ~ /^@/) {
    print
}
