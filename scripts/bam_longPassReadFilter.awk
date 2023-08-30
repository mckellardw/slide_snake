# filter reads in a bam file by length (short pass filter)

# Usage:
## samtools view -h input.bam | awk -v min_length=42 -f bam_longPassReadFilter.awk > output.sam

BEGIN { 
    FS = OFS = "\t" 
}

(length($10) > min_length || $1 ~ /^@/) {
    print
}
