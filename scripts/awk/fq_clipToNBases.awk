# fq_clipToNBases.awk
##  Clip reads to be under a specified length

# Usage:
## awk -v maxLength=650 -f fq_clipToNBases.awk input.fastq > clipped.fastq

NR % 4 == 2 {
    if (length($0) > maxLength) {
        # Clip the sequence to the first maxLength bases
        seq = substr($0, 1, maxLength)
        print seq
    } else {
        print $0
    }
}

NR % 4 == 0 {
    if (length(seq) > 0) {
        # Clip the quality string to match the sequence length
        qual = substr($0, 1, length(seq))
        print qual
    } else {
        print $0
    }
}

NR % 4 != 2 && NR % 4 != 0 {
    print $0
}