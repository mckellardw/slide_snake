# bam_chr2tag.awk
## Add the chromosome to a tag (deafult: "GN)

# Usage:
## samtools view -h input.bam | awk -v tag=GN -f bam_removeTags.awk > output.sam


BEGIN {
    OFS = "\t"
}

# Keep header
/^@/ {
    print
    next
}

{
    if ($3 != "*" && $3 != "=")
        $12 = $12 "\t"tag":Z:" $3
    else
        $12 = $12 "\t"tag":Z:NA"
    print
}
