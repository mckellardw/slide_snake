# bam_lowPassMAPQFilter.awk
## Filter reads in a BAM file based on MAPQ values (keeps reads with MAPQ < threshold)

# Usage:
## samtools view -h input.bam | awk -v quality=30 -f bam_lowPassMAPQFilter.awk > output.sam

BEGIN {
    FS = OFS = "\t"
}

# Keep header lines and reads with low MAPQ
{
    if ($1 ~ /^@/ || $5 < quality) {
        print
    }
}
