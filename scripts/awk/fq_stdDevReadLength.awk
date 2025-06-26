# fq_stdDevReadLength.awk
## Calculate the standard deviation of read lengths in a FASTQ file

# Usage:
## zcat input.fq.gz | awk -f fq_stdDevReadLength.awk

BEGIN {
    sum = 0
    sum_sq = 0
    count = 0
    mean = 0
}

# Process sequence lines (every 4th line starting from line 2)
{
    if (NR % 4 == 2) {
        len = length($0)
        sum += len
        sum_sq += len * len
        count++
    }
}

END {
    if (count > 0) {
        mean = sum / count
        variance = sum_sq / count - mean * mean
        std_dev = sqrt(variance)
        print std_dev
    } else {
        print "No reads found."
    }
}
