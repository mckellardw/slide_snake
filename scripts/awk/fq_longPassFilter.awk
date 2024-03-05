# fq_longPassFilter.awk
##  Filter out reads that are *longer* than the specified length (`minLength`)

# Usage:
## awk -v minLength=650 -f fq_longPassFilter.awk input.fastq > longReads.fastq


{
    pos = NR % 4;
    if (pos == 2) {
        if (length($0) <= minLength) {
            print $0;
        }
    } else {
        print $0;
    }
}