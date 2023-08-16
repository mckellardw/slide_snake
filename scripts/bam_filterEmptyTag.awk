# bam_filterEmptyTag.awk
##  Remove reads/entries which are missing a tag (useful for reads which don't have cell barcodes or UIs associated)

# filter_tag.awk

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
