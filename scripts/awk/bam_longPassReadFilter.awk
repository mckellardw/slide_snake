# bam_longPassReadFilter.awk
## Filter reads in a BAM file by length (long pass filter - keeps reads longer than threshold)

# Usage:
## samtools view -h input.bam | awk -v min_length=42 -f bam_longPassReadFilter.awk > output.sam

BEGIN { 
    FS = OFS = "\t" 
}

# Keep header lines and reads longer than min_length
(length($10) > min_length || $1 ~ /^@/) {
    print
}
