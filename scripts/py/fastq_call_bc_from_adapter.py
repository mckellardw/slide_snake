import argparse
import pysam
from Bio import pairwise2
from Bio.Seq import Seq

# Usage:
## python script_name.py --fastq_file fastq_file.fastq --adapters adapter1 adapter2 --whitelist_files whitelist1.txt whitelist2.txt --barcode_positions left right --mismatches 1 1 --error_rate 1


def align_sequences(seq1, seq2, mismatches):
    alignments = pairwise2.align.globalxx(
        seq1, seq2, score_only=False, one_alignment_only=True
    )
    for a in alignments:
        if a.score >= (len(seq2) - mismatches):
            return a
    return None


def extract_barcodes(
    fastq_file, adapters, whitelist_files, barcode_positions, mismatches, error_rate
):
    barcode_dict = {f: [] for f in whitelist_files}

    with pysam.FastxFile(fastq_file) as fastq:
        for read in fastq:
            for adapter, whitelist_file, position, mismatches in zip(
                adapters, whitelist_files, barcode_positions, mismatches
            ):
                with open(whitelist_file, "r") as whitelist:
                    whitelist = set(whitelist.read().splitlines())
                    if position == "left":
                        barcode_seq = read.sequence[: len(adapter)]
                    else:
                        barcode_seq = read.sequence[-len(adapter) :]
                    alignment = align_sequences(
                        Seq(barcode_seq), Seq(adapter), mismatches
                    )
                    if alignment:
                        barcode = alignment.seqA.tostring()
                        if barcode in whitelist:
                            barcode_dict[whitelist_file].append(barcode)
                            break

    with open("output.tsv", "w") as outfile:
        for whitelist_file, barcodes in barcode_dict.items():
            outfile.write("\t".join(barcodes) + "\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Extract barcodes from FASTQ file.")
    parser.add_argument("--fastq_file", required=True, help="Path to the FASTQ file.")
    parser.add_argument(
        "--adapters", nargs="+", required=True, help="List of adapter sequences."
    )
    parser.add_argument(
        "--whitelist_files", nargs="+", required=True, help="List of whitelist files."
    )
    parser.add_argument(
        "--barcode_positions",
        nargs="+",
        required=True,
        help="List of barcode positions (left or right).",
    )
    parser.add_argument(
        "--mismatches",
        nargs="+",
        type=int,
        required=True,
        help="Number of mismatches allowed in the adapter sequence.",
    )
    parser.add_argument(
        "--error_rate",
        type=int,
        default=1,
        help="Error rate for barcode matching (default: 1 error per 10 bases).",
    )
    args = parser.parse_args()

    extract_barcodes(
        args.fastq_file,
        args.adapters,
        args.whitelist_files,
        args.barcode_positions,
        args.mismatches,
        args.error_rate,
    )
