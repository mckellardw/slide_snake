#!/bin/bash

# Function to print usage information
print_usage() {
    echo "Usage: $0 [options]"
    echo "Options:"
    echo "  -i, --input     Input BAM file (required)"
    echo "  -o, --output    Output BAM file (required)"
    echo "  -t, --tag       Tag label to use (required, must be 2 characters)"
    echo "  -h, --help      Show this help message"
    echo ""
    echo "Example:"
    echo "  $0 -i input.bam -o output.bam -t CR"
    echo "  $0 --input input.bam --output output.bam --tag CR"
}

# Initialize variables
input_bam=""
output_bam=""
tag_label=""

# Parse command line arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        -i|--input)
            input_bam="$2"
            shift 2
            ;;
        -o|--output)
            output_bam="$2"
            shift 2
            ;;
        -t|--tag)
            tag_label="$2"
            shift 2
            ;;
        -h|--help)
            print_usage
            exit 0
            ;;
        *)
            echo "Error: Unknown parameter $1"
            print_usage
            exit 1
            ;;
    esac
done

# Validate required parameters
if [ -z "$input_bam" ] || [ -z "$output_bam" ] || [ -z "$tag_label" ]; then
    echo "Error: Missing required parameters"
    print_usage
    exit 1
fi

# Validate input file exists
if [ ! -f "$input_bam" ]; then
    echo "Error: Input BAM file does not exist: $input_bam"
    exit 1
fi

# Validate tag label format (must be two characters)
if [ ${#tag_label} -ne 2 ]; then
    echo "Error: Tag label must be exactly two characters (e.g., CR, XC)"
    exit 1
fi

# Convert BAM to SAM, add chromosome tags, and convert back to BAM
samtools view -h "$input_bam" \
| awk -v tag="$tag_label" '
BEGIN {OFS="\t"}
{
    if ($0 ~ /^@/) {
        # Print header lines unchanged
        print $0
    } else {
        # For alignment lines, add chromosome as a tag
        # Get chromosome from the third field of the SAM format
        chr=$3
        # Add tag with chromosome information
        print $0 "\t" tag ":Z:" chr
    }
}' \
| samtools view -bS - \
> "$output_bam"
