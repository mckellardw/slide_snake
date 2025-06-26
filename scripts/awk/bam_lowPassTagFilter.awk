# bam_lowPassTagFilter.awk
## Filter reads in a BAM file based on tag values (keeps reads with tag values <= threshold)

# Usage:
## samtools view -h input.bam | awk -v tag=NH -v threshold=1 -f bam_lowPassTagFilter.awk > output.sam

BEGIN {
    FS = OFS = "\t"
}

# Skip header lines
/^@/ { 
    print
    next 
}

# Filter reads based on tag values
{
    for (i = 12; i <= NF; i++) {
        if ($i ~ tag":i:") {
            split($i, arr, ":")
            if (arr[3] <= threshold) {
                print
            }
            break
        }
    }
}
