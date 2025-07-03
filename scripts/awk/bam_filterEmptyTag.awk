# bam_filterEmptyTag.awk
## Remove reads/entries which are missing a tag (useful for reads which don't have cell barcodes or UMIs associated)

# Usage:
## samtools view -h input.bam | awk -v tag=CB -f bam_filterEmptyTag.awk > output.sam

BEGIN {
    OFS = "\t"
}

# Keep header lines
/^@/ {
    print
    next
}

# Filter out reads with empty or missing tags
{
    found = 0
    for(i = 1; i <= NF; i++) {
        if ($i ~ "^"tag":Z:-$") {
            nextfile
        }
        if ($i ~ "^"tag":Z:") {
            found = 1
        }
    }
    if (found == 0) {
        nextfile
    }
    print
}
