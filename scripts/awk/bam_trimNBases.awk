# bam_trimNBases.awk
## Trim N bases from either side of a read sequence in a BAM file

# Usage:
## samtools view -h input.bam | awk -v N=5 -v option=start -f bam_trimNBases.awk > output.sam
## Options: start, end, both

BEGIN { 
    FS = OFS = "\t" 
}

# Process alignment records (skip header lines)
{
    if ($0 !~ /^@/) {
        if (option == "start") {
            $10 = substr($10, N + 1)        # Trim sequence from start
            $11 = substr($11, N + 1)        # Trim quality from start
        } else if (option == "end") {
            $10 = substr($10, 1, length($10) - N)    # Trim sequence from end
            $11 = substr($11, 1, length($11) - N)    # Trim quality from end
        } else if (option == "both") {
            $10 = substr($10, N + 1, length($10) - 2*N)  # Trim from both ends
            $11 = substr($11, N + 1, length($11) - 2*N)  # Trim from both ends
        }
    }
    print $0
}
