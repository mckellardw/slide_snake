# bam_filterMissingTag.awk
## Remove reads/entries that are completely missing a specified tag

# Usage:
## samtools view -h input.bam | awk -v tag=CB -f bam_filterMissingTag.awk > output.sam

BEGIN {
    OFS = "\t"
}

# Keep header lines
/^@/ {
    print
    next
}

# Keep reads that have the specified tag
{
    tag_found = 0
    for (i = 1; i <= NF; i++) {
        if ($i ~ "^"tag":") {
            tag_found = 1
            break
        }
    }
    if (tag_found) {
        print
    } else {
        next
    }
}