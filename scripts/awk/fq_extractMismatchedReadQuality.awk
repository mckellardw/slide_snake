# fq_extractMismatchedReadQuality.awk
##  Extract reads which have mismatched read and quality string lengths

# Usage:
## awk -f fq_extractMismatchedReadQuality.awk input.fastq > mismatchedReads.fastq

{
    pos = NR % 4;
    if (pos == 1) {
        h = $0;
    } else if (pos == 2) {
        s = $0;
    } else if (pos == 3) {
        c = $0;
    } else if (pos == 0) {
        q = $0;
        if (length(q) != length(s)) {
            printf("%s\n%s\n%s\n%s\n", h, s, c, q);
        }
    }
}
