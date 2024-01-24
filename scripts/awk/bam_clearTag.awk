# bam_clearTag.awk
##  Removes a tag from reads/entries

# Usage:
## samtools view -h input.bam | awk -v tag=AS -f bam_removeTags.awk > output.sam

BEGIN {
    OFS = "\t"
}

/^@/ {
    print
    next
}

{
    output = ""
    for(i = 1; i <= NF; i++) {
        if ($i !~ "^"tag":") {
            output = (output == "" ? $i : output OFS $i)
        }
    }
    print output
}