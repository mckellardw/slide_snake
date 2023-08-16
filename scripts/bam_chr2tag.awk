# bam_chr2tag.awk
## Add the chromosome to the "BT" tag

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
        $12 = $12 "\tBT:Z:" $3
    else
        $12 = $12 "\tBT:Z:NA"
    print
}
