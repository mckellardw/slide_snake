# bam_extractTag.awk
##  Write entries for a read ID and a tag to a tab-delimited output file

# Usage:
## samtools view -h input.bam | awk -v tag=CB -f bam_extractTag.awk > output.tsv

BEGIN {
    OFS = "\t"
}

/^@/ {
    read_id = $0   # Save the read ID
    next
}

{
    # Loop through each column to find the desired tag
    for (i = 1; i <= NF; i++) {
        split($i, field, ":")  # Split the column by ':' delimiter
        if (field[1] == tag) {  
            print read_id, field[3]  # Print read ID and tag contents
            break
        }
    }
}
