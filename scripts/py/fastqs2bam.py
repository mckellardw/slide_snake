import argparse
import pysam
import gzip
import tempfile
import os

# Usage:
## python scripts/py/fastqs2bam.py read1.fastq.gz read2.fastq.gz output.bam --r1_seq_tag1 R1 --r1_qual_tag Q1 --r2_qual_tag Q2

# Ref:
## https://readthedocs.org/projects/pysam/downloads/pdf/latest/
## 1.3.4 Creating BAM/CRAM/SAM files from scratch


def fastq_to_bam(r1_fq, r2_fq, out_bam, r1_seq_tag="R1", r1_qual_tag="Q1"):
    """
    Converts paired FASTQ files to an unaligned BAM file, adding barcode information from Read 1 as a tag and quality strings as tags.

    Parameters
    ----------
    r1_fq : str
        Path to the FASTQ file containing Read 1 sequences.
    r2_fq : str
        Path to the FASTQ file containing Read 2 sequences.
    out_bam : str
        Path to the output BAM file.
    r1_seq_tag : str, optional
        Tag name for the barcode information from Read 1, default is 'R1'.
    r1_qual_tag : str, optional
        Tag name for the quality string of Read 1, default is 'Q1'.

    Returns
    -------
    None
    """
    # Create a new BAM file
    header = {"HD": {"VN": "1.0"}, "SQ": []}  # Minimal header
    with pysam.AlignmentFile(out_bam, "wb", header=header) as outfile:
        # Check if the files are gzipped and decompress if necessary
        if r1_fq.endswith(".gz"):
            with gzip.open(r1_fq, "rt") as f:
                r1_fq_content = f.read()
            with tempfile.NamedTemporaryFile(delete=False) as temp_read1:
                temp_read1.write(r1_fq_content.encode())
                r1_fq = temp_read1.name
        if r2_fq.endswith(".gz"):
            with gzip.open(r2_fq, "rt") as f:
                r2_fq_content = f.read()
            with tempfile.NamedTemporaryFile(delete=False) as temp_read2:
                temp_read2.write(r2_fq_content.encode())
                r2_fq = temp_read2.name

        # Open FASTQ files
        with pysam.FastxFile(r1_fq) as read1, pysam.FastxFile(r2_fq) as read2:
            for r1, r2 in zip(read1, read2):
                # Create a new AlignedSegment object
                aln = pysam.AlignedSegment()
                aln.query_name = r1.name
                aln.query_sequence = r2.sequence
                # Convert the quality string to a byte array and assign it
                aln.query_qualities = pysam.qualitystring_to_array(r2.quality)

                # Add barcode information as a tag
                aln.set_tag(r1_seq_tag, r1.sequence)

                # Add quality strings as tags
                aln.set_tag(r1_qual_tag, r1.quality)

                # Write to BAM file
                outfile.write(aln)

        # Clean up temporary files
        if r1_fq.endswith(".gz"):
            os.remove(r1_fq)
        if r2_fq.endswith(".gz"):
            os.remove(r2_fq)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Convert paired FASTQ files to an unaligned BAM file."
    )
    parser.add_argument(
        "r1_fq", help="Path to the FASTQ file containing Read 1 sequences."
    )
    parser.add_argument(
        "r2_fq", help="Path to the FASTQ file containing Read 2 sequences."
    )
    parser.add_argument("out_bam", help="Path to the output BAM file.")
    parser.add_argument(
        "--r1_seq_tag",
        default="R1",
        help="Tag name for the barcode information from Read 1, default is 'R1'.",
    )
    parser.add_argument(
        "--r1_qual_tag",
        default="Q1",
        help="Tag name for the quality string of Read 1, default is 'Q1'.",
    )
    args = parser.parse_args()
    fastq_to_bam(
        args.r1_fq,
        args.r2_fq,
        args.out_bam,
        args.r1_seq_tag,
        args.r1_qual_tag,
    )
