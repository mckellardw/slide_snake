import argparse
import pysam

# Usage:
## python tsv2tag.py path/to/in_bam path/to/tsv_file path/to/out_bam --tag GN


def parse_tsv(tsv_file):
    """Parse the TSV file and return a dictionary mapping read IDs to gene assignments."""
    read_to_tag = {}
    with open(tsv_file, "r") as file:
        for line in file:
            # TODO: abstract this out...
            read_id, status, n_targets, gene = line.strip().split("\t")
            read_to_tag[read_id] = gene
    return read_to_tag


def add_gene_tag_to_bam(in_bam, tsv_file, out_bam, tag="GN"):
    """Add gene assignment from TSV to BAM as a tag and write to a new BAM file."""
    read_to_tag = parse_tsv(tsv_file)
    with pysam.AlignmentFile(in_bam, "rb") as bam_in, pysam.AlignmentFile(
        out_bam, "wb", template=bam_in
    ) as bam_out:
        for read in bam_in:
            if read.query_name in read_to_tag:
                read.set_tag(tag, read_to_tag[read.query_name])
            bam_out.write(read)


def main():
    parser = argparse.ArgumentParser(
        description="Add gene assignment from TSV to BAM as a tag."
    )

    parser.add_argument("in_bam", help="Path to the BAM file.")
    parser.add_argument("tsv_file", help="Path to the TSV file.")
    parser.add_argument("out_bam", help="Path to the output BAM file.")
    parser.add_argument(
        "--tag", default="GN", help='Tag to use for gene assignment. Default is "GN".'
    )
    args = parser.parse_args()

    add_gene_tag_to_bam(args.in_bam, args.tsv_file, args.out_bam, args.tag)


if __name__ == "__main__":
    main()
