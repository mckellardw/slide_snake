# fq_shortPassFilter.awk
##  Filter out reads that are *shorter* than the specified length (`maxLength`)

# Usage:
## awk -v maxLength=650 -f fq_shortPassFilter.awk input.fastq > shortReads.fastq


{
    pos = NR % 4;
    if (pos == 2) {
        if (length($0) >= maxLength) {
            print $0;
        }
    } else {
        print $0;
    }
}