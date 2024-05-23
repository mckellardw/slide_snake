import argparse
import pysam
from Bio import Align
from Bio.Seq import Seq

# Usage:
# SlideSeq
""" 
python scripts/py/fastq_call_bc_from_adapter.py \
    --fq_in data/test_seeker/heart_ctrl_ont.fq.gz \
    --tsv_out barcodes.tsv \
    --adapters ACTGACTG TGCATGCAT \
    --barcode_lengths 10 10 \
    --barcode_positions left right \
    --mismatches 1 1 \
    --error_rate 1 1
"""

# Visium
""" 
python scripts/py/fastq_call_bc_from_adapter.py \
    --fq_in data/test_visium/vy3C_1k_ont.fq.gz \
    --tsv_out sandbox/barcodes.tsv \
    --adapters CTACACGACGCTCTTCCGATCT \
    --barcode_lengths 16 \
    --barcode_positions right \
    --mismatches 2 
"""


def align_sequences(seq, adapter_seq, aligner, mismatches):
    alignments = aligner.align(seq, adapter_seq)

    # Get best alignment
    alignment = alignments[0] if alignments else None

    return alignment
    # for a in alignments:
    #     print(a.score)
    #     if a.score >= (len(adapter_seq) - mismatches):
    #         return a
    # return None


def extract_barcodes(
    fq_in, tsv_out, adapters, barcode_positions, barcode_lengths, mismatches
):
    bc_match_count = 0
    bc_missing_count = 0
    no_bc_count = 0
    read_count = 0

    # Set up aligner
    # TODO- these params need to be better optimized...
    aligner = Align.PairwiseAligner()
    aligner.mode = "local"  # Use 'local' for local alignment
    aligner.match_score = 4  # Match score
    aligner.mismatch_score = -0.5  # Mismatch score
    aligner.open_gap_score = -6  # Gap opening penalty
    aligner.extend_gap_score = -6  # Gap extension penalty

    L_RANGES = list(zip(adapters, barcode_positions, barcode_lengths, mismatches))

    with pysam.FastxFile(fq_in) as fastq:
        with open(tsv_out, "w") as outfile:

            for read in fastq:
                if read_count % 1000000 == 0:
                    print(f"{read_count} reads processed...")
                bc_to_write = []
                for adapter, position, length, mismatches in L_RANGES:
                    alignment = align_sequences(
                        Seq(read.sequence), Seq(adapter), aligner, mismatches
                    )
                    # if alignment.score < 45: #debug
                    #     print(alignment)
                    #     print(alignment.score)

                    if alignment.score > 2 * len(adapter):
                        start = alignment.aligned[0][0][0]
                        end = alignment.aligned[0][0][1]
                        # print(start)
                        # print(end)
                        if position == "left":
                            barcode_seq = read.sequence[(start - length) : start]
                        elif position == "right":
                            barcode_seq = read.sequence[end : (end + length)]
                        else:
                            print(f"Incorrect barcode position [{position}] specified!")

                        # print(barcode_seq)
                        bc_to_write.append(barcode_seq)
                    else:
                        bc_to_write.append("")
                # end alignment

                # Write barcode to file
                if bc_to_write.count("") == 0:
                    bc_match_count += 1
                    # outfile.write(f"{read.name}\t{'\t'.join(bc_to_write)}\n")
                    outfile.write("{}\t{}\n".format(read.name, "\t".join(bc_to_write)))
                elif bc_to_write.count("") == len(adapters):
                    no_bc_count += 1
                else:
                    bc_missing_count += 1
                read_count += 1
            # end fq iterator
        # end tsv open
    # end fq open
    return bc_match_count, bc_missing_count, no_bc_count


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Extract barcodes from FASTQ file and write them in the header."
    )
    parser.add_argument("--fq_in", required=True, help="Path to the input FASTQ file.")
    parser.add_argument(
        "--tsv_out", required=True, help="Path to the output FASTQ file."
    )
    parser.add_argument(
        "--adapters", nargs="+", required=True, help="List of adapter sequences."
    )
    parser.add_argument(
        "--barcode_positions",
        nargs="+",
        required=True,
        help="List of barcode positions (left or right).",
    )
    parser.add_argument(
        "--barcode_lengths",
        nargs="+",
        type=int,
        required=True,
        help="List of barcode lengths.",
    )
    parser.add_argument(
        "--mismatches",
        nargs="+",
        type=int,
        required=True,
        help="Number of mismatches allowed in the adapter sequence.",
    )
    args = parser.parse_args()

    print(
        f"input fastq:        {args.fq_in}\n"
        f"output tsv:         {args.tsv_out}\n"
        f"adapter sequences:  {args.adapters}\n"
        f"barcode positions:  {args.barcode_positions}\n"
        f"barcode lengths:    {args.barcode_lengths}\n"
        f"mismatches allowed: {args.mismatches}\n"
    )
    print("Running...")
    print("")
    bc_match_count, bc_missing_count, no_bc_count = extract_barcodes(
        args.fq_in,
        args.tsv_out,
        args.adapters,
        args.barcode_positions,
        args.barcode_lengths,
        args.mismatches,
    )
    print("")
    print(
        f"# Reads w/ all barcodes found:            {bc_match_count}\n"
        f"# Reads w/ at least one barcode missing:  {bc_missing_count}\n"
        f"# Reads w/ no barcodes                    {no_bc_count}\n"
    )
