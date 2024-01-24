# Remove alignment info from a .bam, but keep the tags
BEGIN { FS = OFS = "\t" }

!/^@/ {
    $3 = "*"
    $4 = "0"
    $5 = "0"
    $6 = "*"
    $7 = "*"
    $8 = "0"
    $9 = "0"
}

{
    print
}