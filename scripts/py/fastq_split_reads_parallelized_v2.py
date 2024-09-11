# fastq_split_reads_parallelized_v2.py

import argparse
import pysam
import os
import time

# import gzip
from Bio import Align

# import string


# Usage:
# SlideSeq
""" 

"""

# Visium
""" 
python scripts/py/fastq_split_reads_parallelized_v2.py \

"""

# microST
"""

"""

# Create a translation table mapping 'ACTG' to 'TGAC'
tab = str.maketrans("ACTG", "TGAC")


def currentTime():
    return time.strftime("%D | %H:%M:%S", time.localtime())


def parse_args():
    parser = argparse.ArgumentParser(
        description="Splits reads in a FASTQ file at a specified sequence, allowing for mismatches."
    )
    parser.add_argument("--fq_in", type=str, help="The path to the FASTQ file.")
    parser.add_argument(
        "--anchor_seq",
        type=str,
        default="CTACACGACGCTCTTCCGATCT",
        help="The primary sequence to search for in the reads.",
    )
    parser.add_argument(
        "--split_seq",
        type=str,
        default=8 * "T",
        help="A harder to find sequence on which you would like to chop the reads (typically poly(T))",
    )
    parser.add_argument(
        "--split_offset",
        type=int,
        default=28,
        help="Positional offset from the 3-prime end of the `split_seq` alignment from where the read should be split. Positive values split toward the 3-prime end, negative toward the 5-prime",
    )
    parser.add_argument(
        "--max_offset",
        type=int,
        default=200,
        help="Maximum scan distance from the 3' end of the anchor alignment",
    )
    parser.add_argument(
        "--max_errors",
        type=int,
        default=2,
        help="The maximum allowed error rate for the anchor sequence matching.",
    )
    parser.add_argument(
        "--threads", type=int, default=1, help="Number of threads to use."
    )

    args = parser.parse_args()

    return args


def reverse_complement(seq):
    # Use the translation table to find the complement of each base
    # and then reverse the sequence
    return seq.translate(tab)[::-1]


def find_and_split_reads(
    fq_in,
    anchor_seq,
    split_seq,
    split_offset,
    max_offset,
    max_errors=2,
    min_read_length=12,
):
    """
    Finds a defined sequence [split_seq] in the middle of reads from a FASTQ file, allowing for mismatches,
    splits the reads at the identified sequence, and writes each split read into separate FASTQ files.

    Parameters
    ----------
    fq_in : str
        The path to the FASTQ file.
    split_seq : str
        The sequence to search for in the reads.
    max_errors : float
        The maximum allowed error rate for the sequence matching.
    """

    output_prefix = fq_in.strip(".fq.gz")

    read_counter = 0
    too_short_counter = 0  # filtered b/c too short counter
    unaligned_counter = 0
    # update_count=1000000 # how often to print updates [for debugging]

    aligner = Align.PairwiseAligner()
    aligner.mode = "local"  # Use 'local' for local alignment
    aligner.match_score = 4  # Match score
    aligner.mismatch_score = -0.5  # Mismatch score
    aligner.open_gap_score = -6  # Gap opening penalty
    aligner.extend_gap_score = -6  # Gap extension penalty

    with pysam.FastxFile(fq_in) as input_fastq, open(
        f"{output_prefix}_R1.fq", mode="w"
    ) as output_fastq1, open(f"{output_prefix}_R2.fq", mode="w") as output_fastq2, open(
        f"{output_prefix}_ambiguous.fq", mode="w"
    ) as output_ambiguous:

        for read in input_fastq:
            # Align anchor sequence
            anchor_alignments = aligner.align(read.sequence, anchor_seq)
            if not anchor_alignments:
                # write to ambiguous.fq
                unaligned_counter += 1
                output_ambiguous.write(
                    f"@{read.name}\n{read.sequence}\n+\n{read.quality}\n"
                )
                continue

            best_anchor_alignment = anchor_alignments[0]
            anchor_end = best_anchor_alignment.aligned[0][0][1]

            # Scan for split sequence
            scan_start = min(0, anchor_end)
            scan_end = max(len(read.sequence), anchor_end + max_offset)
            scan_region = read.sequence[scan_start:scan_end]

            split_alignments = aligner.align(scan_region, split_seq)
            if not split_alignments:
                # write to ambiguous.fq
                unaligned_counter += 1
                output_ambiguous.write(
                    f"@{read.name}\n{read.sequence}\n+\n{read.quality}\n"
                )
                continue

            best_split_alignment = split_alignments[0]
            split_start = scan_start + best_split_alignment.aligned[0][0][0]
            split_end = scan_start + best_split_alignment.aligned[0][0][1]

            # Apply split offset
            # split_point = split_end + split_offset
            split_point = split_end

            # Split the read
            part1_sequence = read.sequence[:split_point]
            part1_quality = read.quality[:split_point]
            part2_sequence = read.sequence[split_point:]
            part2_quality = read.quality[split_point:]

            if (
                len(part1_sequence) > min_read_length
                and len(part2_sequence) > min_read_length
            ):
                output_fastq1.write(
                    f"@{read.name}\n{part1_sequence}\n+\n{part1_quality}\n"
                )
                output_fastq2.write(
                    f"@{read.name}\n{reverse_complement(part2_sequence)}\n+\n{part2_quality[::-1]}\n"
                )
                read_counter += 1
            else:
                # write to ambiguous.fq
                output_ambiguous.write(
                    f"@{read.name}\n{read.sequence}\n+\n{read.quality}\n"
                )
                too_short_counter += 1

            # read count update [for debugging]
            # if read_counter % update_count == 0:
            #     print(f"{currentTime()} - {read_counter} reads processed")

        print(
            f"{read_counter} reads written to {output_prefix}_R1.fq and {output_prefix}_R2.fq.\n"
            f"{too_short_counter} reads removed (shorter than {min_read_length} bases).\n"
            f"{unaligned_counter} reads removed (anchor/split sequences not found).\n"
        )


if __name__ == "__main__":
    args = parse_args()
    # Print run settings for log files ----
    print(
        f"input fastq:          {args.fq_in}\n"
        f"anchor sequence:      {args.anchor_seq}\n"
        f"split sequence:       {args.split_seq}\n"
        f"split offset:         {args.split_offset}\n"
        f"max offset:           {args.max_offset}\n"
        f"max anchor errors:    {args.max_errors}\n"
        f"threads:              {args.threads}\n"
    )
    print(f"{currentTime()} - Running...")
    print("")

    if args.threads == 1:
        find_and_split_reads(
            fq_in=args.fq_in,
            anchor_seq=args.anchor_seq,
            split_seq=args.split_seq,
            split_offset=args.split_offset,
            max_offset=args.max_offset,
            max_errors=args.max_errors,
        )
    elif args.threads > 1:
        # Source: https://superfastpython.com/multiprocessing-pool-for-loop/
        import multiprocessing

        # Split .fq file into {n_core} chunks
        os.system(
            # Custom script to split .fq w/ sed in parallel
            f""" 
            python scripts/py/splitNfqs.py {args.fq_in} {args.threads} {args.threads}
            """
        )

        # Trim each chunked file, in parallel
        chunked_fqs_in = [
            f"{args.fq_in.strip('.fq.gz')}_{str(n).zfill(3)}.fq"
            for n in list(range(1, args.threads + 1))
        ]

        chunked_fqs_out_R1 = [
            f"{args.fq_in.strip('.fq.gz')}_{str(n).zfill(3)}_R1.fq"
            for n in list(range(1, args.threads + 1))
        ]
        chunked_fqs_out_R2 = [
            f"{args.fq_in.strip('.fq.gz')}_{str(n).zfill(3)}_R2.fq"
            for n in list(range(1, args.threads + 1))
        ]
        chunked_fqs_out_ambiguous = [
            f"{args.fq_in.strip('.fq.gz')}_{str(n).zfill(3)}_ambiguous.fq"
            for n in list(range(1, args.threads + 1))
        ]

        # multiprocessing
        items = [
            (
                f"{args.fq_in.strip('.fq.gz')}_{str(n).zfill(3)}.fq",
                args.anchor_seq,
                args.split_seq,
                args.split_offset,
                args.max_offset,
                args.max_errors,
            )
            for n in list(range(1, args.threads + 1))
        ]

        with multiprocessing.Pool(args.threads) as pool:
            multi_out = pool.starmap(find_and_split_reads, items)

        # Merge and compress chunked/trimmed fastqs
        os.system(
            f"""
            cat {' '.join(chunked_fqs_out_R1)} > {args.fq_in.replace('.fq.gz','_R1.fq')}
            cat {' '.join(chunked_fqs_out_R2)} > {args.fq_in.replace('.fq.gz','_R2.fq')}
            cat {' '.join(chunked_fqs_out_ambiguous)} > {args.fq_in.replace('.fq.gz','_ambiguous.fq')}
            """
        )

        if os.path.isfile(args.fq_in.replace(".fq.gz", "_R2.fq.gz")):
            os.system(
                f"""
                rm {' '.join(chunked_fqs_in)} {' '.join(chunked_fqs_out_R1)} {' '.join(chunked_fqs_out_R2)} {' '.join(chunked_fqs_out_ambiguous)}
                """
            )
    else:
        print(f"Value given for '--threads' was not understood. Try again.")
    # end if statement

    # Compress output files
    os.system(
        f"""
        pigz --force -p {args.threads} {args.fq_in.replace('.fq.gz','_R*.fq')}
        pigz --force -p {args.threads} {args.fq_in.replace('.fq.gz','_ambiguous.fq')}
        """
    )
