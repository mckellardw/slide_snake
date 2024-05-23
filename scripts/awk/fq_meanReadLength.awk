# Compute mean read length for a fastq.
## Usage:
#     zcat out/Heart_Control/rRNA/ribodetector/noRibo_R2.fq.gz | head -n 100 | awk -f scripts/awk/fq_meanReadLength.awk 

BEGIN {
    sum = 0;
    count = 0;
}

{
    if (NR % 4 == 2) {
        len = length($0);
        sum += len;
        count++;
    }
}

END {
    mean = sum / count;
    print mean;
}
