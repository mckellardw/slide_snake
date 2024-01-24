# bam_filterEmptyTag.awk
##  Remove reads/entries which are missing a tag (useful for reads which don't have cell barcodes or UIs associated)

# Usage:
## samtools view -h input.bam | awk -v tag=CB -f bam_filterEmptyTag.awk > output.sam

BEGIN {
    OFS = "\t"
}

/^@/ {
    print
    next
}

{
    for(i = 1; i <= NF; i++) {
        if ($i ~ "^"tag":Z:-$") {
            next
        }
    }
    print
}
