import argparse
import pysam
from Bio import Align
from Bio.Seq import Seq
import time

# Usage:
# SlideSeq
""" 
python scripts/py/fastq_call_bc_umi_from_adapter.py \
    --fq_in data/test_seeker/heart_ctrl_ont.fq.gz \
    --tsv_out barcodes.tsv \
    --bc_adapters TCTTCAGCGTTCCCGAGA TCTTCAGCGTTCCCGAGA \
    --bc_lengths 8 6 \
    --bc_offsets 0 \
    --bc_positions left right \
    --bc_mismatches 2 \
    --umi_adapters TCTTCAGCGTTCCCGAGA \
    --umi_lengths 7 \
    --umi_offsets 6 \
    --umi_positions right \
    --umi_mismatches 2 
"""

# Visium
""" 
python scripts/py/fastq_call_bc_umi_from_adapter.py \
    --fq_in data/test_visium/vy3C_1k_ont.fq.gz \
    --tsv_out sandbox/barcodes.tsv \
    --bc_adapters CTACACGACGCTCTTCCGATCT \
    --bc_lengths 16 \
    --bc_offsets 0 \
    --bc_positions right \
    --bc_mismatches 2 \
    --umi_adapters CTACACGACGCTCTTCCGATCT \
    --umi_lengths 12 \
    --umi_offsets 16 \
    --umi_positions right \
    --umi_mismatches 2 
"""


def currentTime():
    return time.strftime("%D | %H:%M:%S", time.localtime())


# Align adapter to read, return alignment
def align_sequences(read, adapter, aligner, mismatches):
    alignments = aligner.align(read, adapter)

    # Get best alignment
    alignment = alignments[0] if alignments else None

    return alignment
    # for a in alignments:
    #     print(a.score)
    #     if a.score >= (len(adapter_seq) - mismatches):
    #         return a
    # return None


# Set up aligner
def init_aligner():
    # TODO- these params need to be better optimized...
    aligner = Align.PairwiseAligner()
    aligner.mode = "local"  # Use 'local' for local alignment
    aligner.match_score = 4  # Match score
    aligner.mismatch_score = -0.5  # Mismatch score
    aligner.open_gap_score = -6  # Gap opening penalty
    aligner.extend_gap_score = -6  # Gap extension penalty

    return aligner


# Simple python version of R rep() function
def rep(val, n):
    return [val for i in range(0, n)]


def main(
    fq_in,
    tsv_out,
    bc_adapters,
    bc_positions,
    bc_lengths,
    bc_offsets,
    bc_mismatches,
    umi_adapters,
    umi_positions,
    umi_lengths,
    umi_offsets,
    umi_mismatches,
):
    bc_match_count = 0
    bc_missing_count = 0
    no_bc_count = 0
    umi_match_count = 0
    umi_missing_count = 0
    no_umi_count = 0
    read_count = 0

    # Set up aligner
    aligner = init_aligner()

    BC_RANGES = list(
        zip(bc_adapters, bc_positions, bc_lengths, bc_offsets, bc_mismatches)
    )
    UMI_RANGES = list(
        zip(umi_adapters, umi_positions, umi_lengths, umi_offsets, umi_mismatches)
    )

    with pysam.FastxFile(fq_in) as fastq:
        with open(tsv_out, "w") as outfile:
            for read in fastq:
                if read_count % 5000000 == 0:
                    print(f"{currentTime()} - {read_count} reads processed...")

                bc_to_write = []
                # TODO- abstract this chunk to a function
                for adapter, position, length, offset, mismatches in BC_RANGES:
                    if len(read.sequence) < len(adapter):
                        continue
                    elif len(read.sequence) < length:
                        continue

                    alignment = align_sequences(
                        Seq(read.sequence), Seq(adapter), aligner, mismatches
                    )
                    # if alignment.score < 45: #debug
                    #     print(alignment)
                    #     print(alignment.score)
                    if alignment is not None:
                        if alignment.score > 2 * len(adapter):
                            start = alignment.aligned[0][0][0]
                            end = alignment.aligned[0][0][1]
                            if position == "left":
                                bc_seq = read.sequence[
                                    (start - offset - length) : (start - offset)
                                ]
                            elif position == "right":
                                bc_seq = read.sequence[
                                    (end + offset) : (end + offset + length)
                                ]
                            else:
                                print(
                                    f"Incorrect barcode position [{position}] specified!"
                                )

                            # Don't write partial barcodes
                            if len(bc_seq) == length:
                                bc_to_write.append(bc_seq)
                            else:
                                bc_to_write.append("-")
                    else:
                        bc_to_write.append("-")
                # end BC alignment

                # Write barcode to file
                ## Only write BC if all positions are found
                if bc_to_write.count("") == 0:
                    bc_match_count += 1

                    # Proceed w/ UMI alignment
                    umi_to_write = []
                    for adapter, position, length, offset, mismatches in UMI_RANGES:
                        alignment = align_sequences(
                            Seq(read.sequence), Seq(adapter), aligner, mismatches
                        )

                        if alignment.score > 2 * len(adapter):
                            start = alignment.aligned[0][0][0]
                            end = alignment.aligned[0][0][1]
                            if position == "left":
                                umi_seq = read.sequence[
                                    (start - offset - length) : (start - offset)
                                ]
                            elif position == "right":
                                umi_seq = read.sequence[
                                    (end + offset) : (end + offset + length)
                                ]
                            else:
                                print(f"Incorrect UMI position [{position}] specified!")

                            # Don't write partial barcodes
                            if len(umi_seq) == length:
                                umi_to_write.append(umi_seq)
                            else:
                                umi_to_write.append("-")
                        else:
                            umi_to_write.append("-")
                    # end UMI alignment

                    outfile.write(
                        "{}\t{}\t{}\n".format(
                            read.name, "\t".join(bc_to_write), "\t".join(umi_to_write)
                        )
                    )
                    # Write UMI to file
                    if umi_to_write.count("") == 0:
                        umi_match_count += 1
                    elif umi_to_write.count("") == len(umi_adapters):
                        no_umi_count += 1
                    else:
                        umi_missing_count += 1
                elif bc_to_write.count("") == len(bc_adapters):
                    no_bc_count += 1
                else:
                    bc_missing_count += 1

                read_count += 1
            # end fq iterator
        # end tsv open
    # end fq open
    return (
        bc_match_count,
        bc_missing_count,
        no_bc_count,
        umi_match_count,
        umi_missing_count,
        no_umi_count,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Extract barcodes from FASTQ file and write them in the header."
    )
    parser.add_argument("--fq_in", required=True, help="Path to the input FASTQ file.")
    parser.add_argument(
        "--tsv_out", required=True, help="Path to the output FASTQ file."
    )
    parser.add_argument(
        "--bc_adapters",
        nargs="+",
        required=True,
        help="List of adapter sequences used for barcode extraction.",
    )
    parser.add_argument(
        "--bc_positions",
        nargs="+",
        required=True,
        help="List of barcode positions (left or right).",
    )
    parser.add_argument(
        "--bc_lengths",
        nargs="+",
        type=int,
        required=True,
        help="List of barcode lengths.",
    )
    parser.add_argument(
        "--bc_offsets",
        nargs="+",
        type=int,
        required=True,
        help="List of offsets from adapter sequence. Keep in mind the `bc_position`",
    )
    parser.add_argument(
        "--bc_mismatches",
        nargs="+",
        type=int,
        required=True,
        help="Number of mismatches allowed in the adapter sequence for BC calling.",
    )
    parser.add_argument(
        "--umi_adapters",
        nargs="+",
        required=True,
        help="List of adapter sequences used for UMI extraction.",
    )
    parser.add_argument(
        "--umi_positions",
        nargs="+",
        required=True,
        help="List of UMI positions (`left` or `right`).",
    )
    parser.add_argument(
        "--umi_lengths",
        nargs="+",
        type=int,
        required=True,
        help="List of UMI lengths.",
    )
    parser.add_argument(
        "--umi_offsets",
        nargs="+",
        type=int,
        required=True,
        help="List of offsets from adapter sequence. Keep in mind the `umi_position`",
    )
    parser.add_argument(
        "--umi_mismatches",
        nargs="+",
        type=int,
        required=True,
        help="Number of mismatches allowed in the adapter sequence for UMI calling.",
    )
    args = parser.parse_args()

    # param checks ------
    if (
        len(args.bc_mismatches) != len(args.bc_adapters)
        and len(args.bc_mismatches) == 1
    ):
        args.bc_mismatches = rep(val=args.bc_mismatches[0], n=len(args.bc_adapters))

    if len(args.bc_offsets) != len(args.bc_adapters) and len(args.bc_offsets) == 1:
        args.bc_offsets = rep(val=args.bc_offsets[0], n=len(args.bc_adapters))

    if (
        len(args.umi_mismatches) != len(args.umi_adapters)
        and len(args.umi_mismatches) == 1
    ):
        args.umi_mismatches = rep(val=args.umi_mismatches[0], n=len(args.umi_adapters))

    if len(args.umi_offsets) != len(args.umi_adapters) and len(args.umi_offsets) == 1:
        args.umi_offsets = rep(val=args.umi_offsets[0], n=len(args.umi_adapters))

    # Print run settings for log files ----
    print(
        f"input fastq:                      {args.fq_in}\n"
        f"output tsv:                       {args.tsv_out}\n"
        f"BC adapter sequence(s):           {args.bc_adapters}\n"
        f"BC position(s):                   {args.bc_positions}\n"
        f"BC offset(s):                     {args.bc_offsets}\n"
        f"BC length(s):                     {args.bc_lengths}\n"
        f"BC adapter mismatches allowed:    {args.bc_mismatches}\n"
        f"UMI adapter sequence(s):          {args.umi_adapters}\n"
        f"UMI position(s):                  {args.umi_positions}\n"
        f"UMI offset(s):                    {args.umi_offsets}\n"
        f"UMI length(s):                    {args.umi_lengths}\n"
        f"UMI adapter mismatches allowed:   {args.umi_mismatches}\n"
    )
    print(f"{currentTime()} - Running...")
    print("")

    # Run main function ------
    (
        bc_match_count,
        bc_missing_count,
        no_bc_count,
        umi_match_count,
        umi_missing_count,
        no_umi_count,
    ) = main(
        args.fq_in,
        args.tsv_out,
        args.bc_adapters,
        args.bc_positions,
        args.bc_lengths,
        args.bc_offsets,
        args.bc_mismatches,
        args.umi_adapters,
        args.umi_positions,
        args.umi_lengths,
        args.umi_offsets,
        args.umi_mismatches,
    )
    print("")
    print(
        f"{currentTime()}\n"
        f"# Reads w/ all barcodes found:            {bc_match_count}\n"
        f"# Reads w/ at least one barcode missing:  {bc_missing_count}\n"
        f"# Reads w/ no barcodes                    {no_bc_count}\n"
        f"# Reads w/ all UMIs found:            {umi_match_count}\n"
        f"# Reads w/ at least one UMI missing:  {umi_missing_count}\n"
        f"# Reads w/ no UMIs                    {no_umi_count}\n"
    )
