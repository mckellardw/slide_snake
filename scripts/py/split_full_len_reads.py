import argparse
import pysam
from Bio import pairwise2
from Bio import Seq
# from Bio.Align import PairwiseAligner
# from Bio.SubsMat import MatrixInfo as matlist
import string

# Create a translation table mapping 'ACTG' to 'TGAC'
tab = str.maketrans("ACTG", "TGAC")

def reverse_complement(seq):
    # Use the translation table to find the complement of each base
    # and then reverse the sequence
    return seq.translate(tab)[::-1]


def find_and_split_reads(fastq_file, output_prefix, sequence, log, max_errors=2):
    """
    Finds an arbitrary sequence in the middle of reads from a FASTQ file, allowing for mismatches,
    splits the reads at the identified sequence, and writes each split read into separate FASTQ files.

    Parameters
    ----------
    fastq_file : str
        The path to the FASTQ file.
    output_prefix : str
        The prefix for the output FASTQ files.
    sequence : str
        The sequence to search for in the reads.
    max_errors : float
        The maximum allowed error rate for the sequence matching.
    """
    read_counter=0

    # Open the FASTQ file
    with pysam.FastxFile(fastq_file) as input_fastq:
        # Open the output FASTQ files
        with open(f"{output_prefix}_R1.fq", mode='w') as output_fastq1, \
             open(f"{output_prefix}_R2.fq", mode='w') as output_fastq2:
            for read in input_fastq:
                # Perform pairwise alignment to find the sequence with allowed mismatches
                alignments = pairwise2.align.localms(
                    read.sequence, 
                    sequence,
                    4, -0.5, -6, -6, #TODO- optimize these...
                    one_alignment_only=True
                )

                # Calculate the error rate for the best alignment
                best_alignment = alignments[0] if alignments else None
                error_rate =  1 - (best_alignment.score / len(sequence)) if best_alignment else float('inf')

                # Check if the error rate is within the allowed limit
                if error_rate <= max_errors:
                    # Split the read into two at the identified sequence
                    sequence_start = best_alignment[3]
                    sequence_end = best_alignment[4]

                    #R1
                    part1_sequence = read.sequence[:sequence_start]
                    part1_quality = read.quality[:sequence_start]

                    #R2
                    part2_sequence = read.sequence[sequence_end:]
                    part2_quality = read.quality[sequence_end:]

                    # Write the split reads to the output FASTQ files
                    output_fastq1.write(
                        # pysam.FastxRecord(name=read.name + "_R1", sequence=part1_sequence, quality=part1_quality)
                        f"@{read.name}\n{part1_sequence}\n+\n{part1_quality}\n"
                    )
                    output_fastq2.write(
                        # pysam.FastxRecord(name=read.name + "_R2", sequence=part2_sequence, quality=part2_quality)
                        f"@{read.name}\n{reverse_complement(part2_sequence)}\n+\n{part2_quality[::-1]}\n" # return rev comp for alignment
                        # f"@{read.name}\n{part2_sequence}\n+\n{part2_quality}"
                    )

                    read_counter+=1
                    #end loop
    
    # Write log
    with open(log, mode='w') as output_log:
        output_log.write(
            f"{read_counter} reads written."
        )

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='Splits reads in a FASTQ file at a specified sequence, allowing for mismatches.'
    )
    parser.add_argument('--fastq_file', 
                        type=str, 
                        help='The path to the FASTQ file.')
    parser.add_argument('--output_prefix', 
                        type=str, 
                        help='The prefix for the output FASTQ files.')
    parser.add_argument('--sequence', 
                        type=str, 
                        help='The sequence to search for in the reads.')
    parser.add_argument('--log', 
                        type=str, 
                        help='The path to the output log file.')
    parser.add_argument('--max_errors', 
                        type=int, 
                        default=2, 
                        help='The maximum allowed error rate for the sequence matching.')
    args = parser.parse_args()

    find_and_split_reads(
        args.fastq_file, 
        args.output_prefix, 
        args.sequence,
        args.log,
        args.max_errors
    )
