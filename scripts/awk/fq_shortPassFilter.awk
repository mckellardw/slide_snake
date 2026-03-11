# fq_shortPassFilter.awk
## Filter out reads that are *shorter* than the specified length (`maxLength`)

# Usage:
## awk -v maxLength=650 -f fq_shortPassFilter.awk input.fastq > shortReads.fastq

# Process FASTQ records (4 lines per record)
{
    header = $0
    getline seq
    getline sep
    getline qual

    if (length(seq) >= maxLength) {
        print header
        print seq
        print sep
        print qual
    }
}
