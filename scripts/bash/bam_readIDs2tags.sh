#!/bin/bash

# Check if the correct number of arguments is provided
if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <input.bam> <output.bam> <tag_prefix>"
    exit 1
fi

input_bam=$1
output_bam=$2
tag_prefix=$3

# Convert BAM to SAM, add tags, and convert back to BAM
samtools view -h $input_bam \
| awk -v prefix=$tag_prefix '
BEGIN {OFS="\t"}
{
    split($1, a, "_");
    if (length(a) >= 2) {
        $1 = a[1];
        print $0 "\t" prefix "CR:" a[length(a)-1] "\t" prefix "UR:" a[length(a)]
    } else {
        print $0
    }
}' \
| samtools view -bS - \
> $output_bam

echo "Tags added to $output_bam"
