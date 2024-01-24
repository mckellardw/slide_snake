# filter reads in a bam file by length (short pass filter)
BEGIN { FS = OFS = "\t" }

(length($10) < max_length || $1 ~ /^@/) {
    print
}
