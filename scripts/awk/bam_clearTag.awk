# bam_clearTag.awk
## Removes a tag from reads/entries

# Usage:
## samtools view -h input.bam | awk -v tag=AS -f bam_clearTag.awk > output.sam

BEGIN {
    OFS = "\t"
}

# Keep header lines
/^@/ {
    print
    next
}

# Remove specified tag from alignment records
{
    output = ""
    for(i = 1; i <= NF; i++) {
        if ($i !~ "^"tag":") {
            output = (output == "" ? $i : output OFS $i)
        }
    }
    print output
}