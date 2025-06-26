# fq_meanReadLength.awk
## Compute mean read length for a FASTQ file

# Usage:
## zcat out/Heart_Control/rRNA/ribodetector/noRibo_R2.fq.gz | head -n 100 | awk -f scripts/awk/fq_meanReadLength.awk

BEGIN {
    sum = 0
    count = 0
}

# Process sequence lines (every 4th line starting from line 2)
{
    if (NR % 4 == 2) {
        len = length($0)
        sum += len
        count++
    }
}

END {
    mean = sum / count
    print mean
}
