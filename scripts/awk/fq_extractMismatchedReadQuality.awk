# fq_extractMismatchedReadQuality.awk
## Extract reads which have mismatched read and quality string lengths

# Usage:
## awk -f fq_extractMismatchedReadQuality.awk input.fastq > mismatchedReads.fastq

# Process FASTQ records and check for length mismatches
{
    pos = NR % 4
    if (pos == 1) {         # Header line
        h = $0
    } else if (pos == 2) {  # Sequence line
        s = $0
    } else if (pos == 3) {  # Separator line
        c = $0
    } else if (pos == 0) {  # Quality line
        q = $0
        if (length(q) != length(s)) {
            printf("%s\n%s\n%s\n%s\n", h, s, c, q)
        }
    }
}
