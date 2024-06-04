import argparse
import pysam

# Usage:
"""
python scripts/py/tsv2tag.py \
    --in_bam sandbox/vy3c.bam \
    --in_tsv barcodes_corrected.tsv \
    --out_bam sandbox/vy3c_cb.bam \
    --readIDColumn 0 \
    --tagColumns 1 \
    --tags CB
"""


# def parse_tsv(in_tsv, readIDcolumn, tagColumn):
#     """Parse the TSV file and return a dictionary mapping read IDs to gene assignments."""
#     read_to_tag = {}
#     with open(in_tsv, "r") as file:
#         for line in file:
#             # TODO: abstract this out...
#             line_as_list = line.strip().split("\t")
#             read_to_tag[readIDcolumn] = line_as_list[tagColumn]
#     return read_to_tag

def parse_tsv(in_tsv, readIDcolumn, tagColumns):
    """Parse the TSV file and return a dictionary mapping read IDs to tag values."""
    read_to_tags = {}
    with open(in_tsv, "r") as file:
        for line in file:
            line_as_list = line.strip().split("\t")
            read_id = line_as_list[readIDcolumn]
            read_to_tags[read_id] = [line_as_list[tagColumn] for tagColumn in tagColumns]
    return read_to_tags


def add_tags_to_bam(in_bam, in_tsv, out_bam, readIDColumn=0, tagColumns=[1], tags=["GN"]):
    """Add tag values from TSV to BAM and write to a new BAM file."""
    read_to_tags = parse_tsv(in_tsv, readIDColumn, tagColumns)
    
    reads_yes_tags=0
    reads_no_tags=0

    with pysam.AlignmentFile(in_bam, "rb", check_sq=False) as bam_in, pysam.AlignmentFile(
        out_bam, "wb", template=bam_in
    ) as bam_out:
        for read in bam_in:
            if read.query_name in read_to_tags:
                reads_yes_tags+=1
                for i, tag_value in enumerate(read_to_tags[read.query_name]):
                    read.set_tag(tags[i], tag_value)
            else:                
                reads_no_tags+=1
            bam_out.write(read)
        return [reads_yes_tags, reads_no_tags]


def main():
    parser = argparse.ArgumentParser(
        description="Add gene assignment from TSV to BAM as a tag."
    )

    parser.add_argument("--in_bam", help="Path to the input BAM file.")
    parser.add_argument("--in_tsv", help="Path to the input TSV file.")
    parser.add_argument("--out_bam", help="Path to the output BAM file.")
    parser.add_argument(
        "--readIDColumn",
        nargs="+",
        type=int,
        default=1,
        help="Which (zero-indexed) column in the TSV file contains the read ID? If passing multiple tag values, separate the list with spaces. Default is 0 (zero-indexed, first column...).",
    )
    parser.add_argument(
        "--tagColumns",
        nargs="+",
        type=int,
        default=1,
        help="Which (zero-indexed) column in the TSV file contains the info that will be added to the output bam? If passing multiple tag values, separate the list with spaces. Default is 1 (zero-indexed, second column...).",
    )
    parser.add_argument(
        "--tags",
        nargs="+",
        default="GN",
        help="Tag(s) to use for gene assignment. If passing multiple tag values, separate the list with spaces. Default is 'GN'.",
    )

    args = parser.parse_args()
    print(
        f"input BAM:                {args.in_bam}\n"
        f"input TSV:                {args.in_tsv}\n"
        f"output BAM:               {args.out_bam}\n"
        f"Read ID column in TSV:    {args.readIDColumn[0]}\n"
        f"Tag column(s) in TSV:     {list(args.tagColumns)}\n"
        f"Tag(s) to add:            {list(args.tags)}\n"
    )
    reads_yes_tags, reads_no_tags = add_tags_to_bam(
        in_bam=args.in_bam,
        in_tsv=args.in_tsv,
        out_bam=args.out_bam,
        readIDColumn=args.readIDColumn[0],
        tagColumns=list(args.tagColumns),
        tags=list(args.tags),
    )
    print("")
    print(
        f"# total reads parsed:         {reads_yes_tags+ reads_no_tags}\n"
        f"# Reads w/ tag(s) found:      {reads_yes_tags}\n"
        f"# Reads w/ no tag(s) found:   {reads_no_tags}\n"
    )



if __name__ == "__main__":
    main()
