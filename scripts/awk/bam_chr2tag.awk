# bam_chr2tag.awk
## Add the chromosome to a tag (default: "GN")

# Usage:
## samtools view -h input.bam | awk -v tag=GN -f bam_chr2tag.awk > output.sam


BEGIN {
    OFS = "\t"
}

# Keep header lines
/^@/ {
    print
    next
}

# Add chromosome information to tag
{
    if ($3 != "*" && $3 != "=")
        $12 = $12 "\t"tag":Z:" $3
    else
        $12 = $12 "\t"tag":Z:NA"
    print
}
