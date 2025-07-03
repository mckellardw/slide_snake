# bam_combineTags.awk
## Combine values from multiple tags into a new tag in a SAM/BAM entry

# Usage:
## samtools view input.bam | awk -v tags=CB,UB -v outtag=XB -f bam_combineTags.awk > output.sam

BEGIN {
    FS = OFS = "\t"
    n_tags = split(tags, tag_list, ",")
}

# Process alignment records (skip header lines)
!/^@/ {
    # Store original fields
    for (i = 1; i <= NF; i++) {
        fields[i] = $i
    }
    nfields = NF

    # Extract values for each tag
    combined = ""
    for (t = 1; t <= n_tags; t++) {
        tagval = ""
        for (i = 12; i <= nfields; i++) {
            split(fields[i], field, ":")
            if (field[1] == tag_list[t]) {
                tagval = field[3]
                break
            }
        }
        combined = combined tagval
    }

    # Print original fields, then append new tag
    for (i = 1; i <= nfields; i++) {
        printf "%s%s", fields[i], (i < nfields ? OFS : "")
    }
    printf "%s%s:Z:%s\n", OFS, outtag, combined
}