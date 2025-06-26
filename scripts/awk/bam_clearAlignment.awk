# bam_clearAlignment.awk
## Remove alignment info from a BAM file, but keep the tags

# Usage:
## samtools view -h input.bam | awk -f bam_clearAlignment.awk > output.sam

BEGIN { 
    FS = OFS = "\t" 
}

# Clear alignment fields for non-header lines
!/^@/ {
    $3 = "*"    # Reference name
    $4 = "0"    # Position
    $5 = "0"    # MAPQ
    $6 = "*"    # CIGAR
    $7 = "*"    # Mate reference name
    $8 = "0"    # Mate position
    $9 = "0"    # Template length
}

# Print all lines (headers and modified alignments)
{
    print
}