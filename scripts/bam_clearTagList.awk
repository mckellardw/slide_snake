# bam_clearTagList.awk
##  Removes a list of tags from reads/entries

# Usage:
## samtools view -h input.bam | awk -v tags="TAG1,TAG2,TAG3" -f bam_removeTags.awk > output.sam

BEGIN {
    OFS = "\t"
    n = split(tags, a, ",")
    for(i = 1; i <= n; i++) {
        taglist[a[i]] = 1
    }
}

/^@/ {
    print
    next
}

{
    output = ""
    for(i = 1; i <= NF; i++) {
        if (!($i in taglist)) {
            output = (output == "" ? $i : output OFS $i)
        }
    }
    print output
}
