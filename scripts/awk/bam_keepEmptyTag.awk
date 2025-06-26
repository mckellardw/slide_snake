# bam_keepEmptyTag.awk
## Keep reads/entries which have an empty tag (useful for filtering reads with missing cell barcodes or UMIs)

# Usage:
## samtools view -h input.bam | awk -v tag=CB -f bam_keepEmptyTag.awk > output.sam

BEGIN {
    OFS = "\t"
}

# Keep header lines
/^@/ {
    print
    next
}

# Keep reads with empty tags
{
    for(i = 1; i <= NF; i++) {
        if ($i ~ "^"tag":Z:-$") {
            print
            next
        }
    }
}
