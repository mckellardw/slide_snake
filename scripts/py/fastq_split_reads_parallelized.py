import argparse
import pysam
import os

# import gzip
# from Bio import pairwise2
# from Bio import Seq
from Bio import Align

# from Bio.Align import PairwiseAligner
# from Bio.SubsMat import MatrixInfo as matlist
import string

# Create a translation table mapping 'ACTG' to 'TGAC'
tab = str.maketrans("ACTG", "TGAC")


def reverse_complement(seq):
    # Use the translation table to find the complement of each base
    # and then reverse the sequence
    return seq.translate(tab)[::-1]


def find_and_split_reads(fq_in, split_seq, log=None, max_errors=2, min_read_length=12):
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

    aligner = Align.PairwiseAligner()
    aligner.mode = "local"  # Use 'local' for local alignment
    aligner.match_score = 4  # Match score
    aligner.mismatch_score = -0.5  # Mismatch score
    aligner.open_gap_score = -6  # Gap opening penalty
    aligner.extend_gap_score = -6  # Gap extension penalty

    # Open the FASTQ file
    with pysam.FastxFile(fq_in) as input_fastq:
        # Open the output FASTQ files
        with open(f"{output_prefix}_R1.fq", mode="w") as output_fastq1, open(
            f"{output_prefix}_R2.fq", mode="w"
        ) as output_fastq2:
            for read in input_fastq:
                # Perform pairwise alignment to find the sequence with allowed mismatches
                alignments = aligner.align(read.sequence, split_seq)

                # Calculate the error rate for the best alignment
                best_alignment = alignments[0] if alignments else None
                error_rate = (
                    1 - (best_alignment.score / len(split_seq))
                    if best_alignment
                    else float("inf")
                )

                # print(best_alignment.aligned[0][0][0])
                # print(best_alignment.aligned[0][0][1])
                # print("")
                # Check if the error rate is within the allowed limit
                if error_rate <= max_errors:
                    # Split the read into two at the identified sequence
                    sequence_start = best_alignment.aligned[0][0][0]
                    sequence_end = best_alignment.aligned[0][0][1]

                    # R1 - include the poly(T) region
                    part1_sequence = read.sequence[:sequence_end]
                    part1_quality = read.quality[:sequence_end]

                    # R2 - include the poly(T) region as well, rev-comp it to a poly(A) tail
                    # part2_sequence = read.sequence[sequence_end:]
                    # part2_quality = read.quality[sequence_end:]
                    part2_sequence = read.sequence[sequence_start:]
                    part2_quality = read.quality[sequence_start:]

                    if (
                        len(part1_sequence) > min_read_length
                        and len(part1_sequence) > min_read_length
                    ):
                        # Write the split reads to the outpu    t FASTQ files
                        output_fastq1.write(
                            # pysam.FastxRecord(name=read.name + "_R1", sequence=part1_sequence, quality=part1_quality)
                            f"@{read.name}\n{part1_sequence}\n+\n{part1_quality}\n"
                        )
                        output_fastq2.write(
                            # pysam.FastxRecord(name=read.name + "_R2", sequence=part2_sequence, quality=part2_quality)
                            f"@{read.name}\n{reverse_complement(part2_sequence)}\n+\n{part2_quality[::-1]}\n"  # return rev comp for alignment
                            # f"@{read.name}\n{part2_sequence}\n+\n{part2_quality}"
                        )
                        read_counter += 1
                    else:
                        too_short_counter += 1
                    # end loop

    if log is not None:
        # Write log
        with open(log, mode="w") as output_log:
            output_log.write(
                f"{read_counter} reads written to {output_prefix}_R1.fq and {output_prefix}_R2.fq.\n"
            )
            output_log.write(
                f"{too_short_counter} reads removed (shorter than {min_read_length} bases).\n"
            )
    else:
        print(
            f"{read_counter} reads written to {output_prefix}_R1.fq and {output_prefix}_R2.fq.\n"
        )
        print(
            f"{too_short_counter} reads removed (shorter than {min_read_length} bases).\n"
        )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Splits reads in a FASTQ file at a specified sequence, allowing for mismatches."
    )

    parser.add_argument("--fq_in", type=str, help="The path to the FASTQ file.")
    parser.add_argument(
        "--split_seq", type=str, help="The sequence to search for in the reads."
    )
    # parser.add_argument('--output_prefix',
    #                     type=str,
    #                     default=None,
    #                     help='The prefix for the output FASTQ files.')
    parser.add_argument(
        "--log", type=str, default="split.log", help="The path to the output log file."
    )
    parser.add_argument(
        "--threads", type=int, default=1, help="Number of threads to use."
    )
    parser.add_argument(
        "--max_errors",
        type=int,
        default=2,
        help="The maximum allowed error rate for the sequence matching.",
    )
    args = parser.parse_args()

    if args.threads == 1:
        find_and_split_reads(
            fq_in=args.fq_in,
            split_seq=args.split_seq,
            # output_prefix = args.output_prefix,
            log=args.log,
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

        # multiprocessing
        items = [
            (
                f"{args.fq_in.strip('.fq.gz')}_{str(n).zfill(3)}.fq",
                args.split_seq,
                # args.output_prefix,
                # args.log,
                args.max_errors,
            )
            for n in list(range(1, args.threads + 1))
        ]

        with multiprocessing.Pool(args.threads) as pool:
            multi_out = pool.starmap(find_and_split_reads, items)

        # TODO- logging

        # Merge and compress chunked/trimmed fastqs
        os.system(
            f"""
            cat {' '.join(chunked_fqs_out_R1)} > {args.fq_in.replace('.fq.gz','_R1.fq')}
            cat {' '.join(chunked_fqs_out_R2)} > {args.fq_in.replace('.fq.gz','_R2.fq')}
            """
        )

        if os.path.isfile(args.fq_in.replace(".fq.gz", "_R2.fq.gz")):
            os.system(
                f"""
                rm {' '.join(chunked_fqs_in)} {' '.join(chunked_fqs_out_R1)} {' '.join(chunked_fqs_out_R2)}
                """
            )
    else:
        print(f"Value given for '--threads' was not understood. Try again.")
    # end if statement

    # Compress output files
    os.system(
        f"""
        pigz --force -p {args.threads} {args.fq_in.replace('.fq.gz','_R*.fq')}
        """
    )
