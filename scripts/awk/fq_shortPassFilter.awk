# fq_shortPassFilter.awk
## Filter out reads that are *shorter* than the specified length (`maxLength`)

# Usage:
## awk -v maxLength=650 -f fq_shortPassFilter.awk input.fastq > shortReads.fastq

# TODO - fix this so that the extra "+" and read IDs are not also written...

# Process FASTQ records
{
    pos = NR % 4
    if (pos == 2) {         # Sequence line
        if (length($0) >= maxLength) {
            print $0
        }
    } else if (pos == 4) {  # Quality line
        if (length($0) >= maxLength) {
            print $0
        }
    } else {                # Header and separator lines
        print $0
    }
}
