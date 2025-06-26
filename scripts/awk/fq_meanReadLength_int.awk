# fq_meanReadLength_int.awk
## Compute mean read length for a FASTQ file (integer output)

# Usage:
## zcat out/Heart_Control/rRNA/ribodetector/noRibo_R2.fq.gz | head -n 100 | awk -f scripts/awk/fq_meanReadLength_int.awk

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
    print int(mean)
}
