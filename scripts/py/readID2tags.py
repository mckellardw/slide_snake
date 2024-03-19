#!/usr/bin/env python3

"""
Script to parse a BAM file, extract the read ID, and move the last two strings in the ID to tags.

Usage:
    python3 readID2tags.py -i input.bam -o output.bam -c CR -y UR

Options:
    -i, --input      Path to the input BAM file.
    -o, --output     Path to the output BAM file.
    -c, --barcode_tag     The tag for the first string in the read ID. Default: CR
    -y, --umi_tag     The tag for the second string in the read ID. Default: CY
"""

import argparse
import pysam

def parse_bam_and_add_tags(input_bam, output_bam, barcode_tag, umi_tag):
    with pysam.AlignmentFile(input_bam, "rb") as infile, \
         pysam.AlignmentFile(output_bam, "wb", template=infile) as outfile:
        
        for read in infile:
            read_id = read.query_name
            parts = read_id.split("_")
            if len(parts) >= 2:
                barcode_value = parts[-2]
                umi_value = parts[-1]
                
                # Extract quality scores
                # qual_scores = read.query_qualities
                
                read.set_tag(barcode_tag, barcode_value)
                read.set_tag(umi_tag, umi_value)
            
            outfile.write(read)

def main():
    parser = argparse.ArgumentParser(
        description="Add tags to BAM file based on read ID."
    )
    parser.add_argument(
        "-i", "--input", 
        required=True, 
        help="Path to the input BAM file."
    )
    parser.add_argument(
        "-o", "--output", 
        required=True, 
        help="Path to the output BAM file."
    )
    parser.add_argument(
        "-c", "--barcode_tag", 
        default="CR", 
        help="The tag for the first string in the read ID [which should contain the cell barcode]. Default: CR"
    )
    parser.add_argument(
        "-y", "--umi_tag", 
        default="UR",
        help="The tag for the second string in the read ID [which should contain the UMI]. Default: UR"
    )
    #TODO:
    # parser.add_argument(
    #     "--barcode_qual_tag", 
    #     default="CY", 
    #     help="Read quality for the barcode. Default: CY"
    # )
    # parser.add_argument(
    #     "--umi_qual_tag", 
    #     default="UY",
    #     help="Read quality for the UI. Default: UY"
    # )
    
    args = parser.parse_args()
    
    parse_bam_and_add_tags(args.input, args.output, args.barcode_tag, args.umi_tag)

if __name__ == "__main__":
    main()
