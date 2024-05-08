import argparse
import os
import pysam
import time
from itertools import groupby
from collections import Counter


def calculate_metrics_bam(bam_file, tags, chunk_size=10000):
    metrics = []
    with pysam.AlignmentFile(bam_file, "rb") as samfile:
        for i, read in enumerate(samfile):
            read_length = read.query_length
            if read_length > 0:
                gc_count = (
                    Counter(read.query_sequence)["G"]
                    + Counter(read.query_sequence)["C"]
                )
                gc_percent = (gc_count / read_length) * 100
                first_base = read.query_sequence[0]
                last_base = read.query_sequence[-1]
                homopolymer_lengths = [
                    len(list(g)) for k, g in groupby(read.query_sequence)
                ]
                longest_homopolymer = max(homopolymer_lengths)
                homopolymer_base = read.query_sequence[
                    homopolymer_lengths.index(longest_homopolymer)
                ]
                alignment_status = {
                    "is_unmapped": read.is_unmapped,
                    "is_proper_pair": read.is_proper_pair,
                    "is_secondary": read.is_secondary,
                }
                # tag_values = {tag: read.get_tag(tag) for tag in tags}
                tag_values = {}
                for tag in tags:
                    try:
                        tag_values[tag] = read.get_tag(tag)
                    except:
                        tag_values[tag] = None

                metrics.append(
                    {
                        "Read_ID": read.query_name,
                        "Read_Length": read_length,
                        "GC_Percent": round(gc_percent, 2),
                        "First_Base": first_base,
                        "Last_Base": last_base,
                        "Longest_Homopolymer": longest_homopolymer,
                        "Homopolymer_Base": homopolymer_base,
                        **alignment_status,
                        **tag_values,
                    }
                )

                if (i + 1) % chunk_size == 0:
                    yield metrics
                    metrics = []
            # fi
        # rof

        if metrics:  # Yield any remaining metrics
            yield metrics


def write_tsv(metrics, tsv_file, force_overwrite=True):
    if not force_overwrite and os.path.exists(tsv_file):
        print(f"File {tsv_file} already exists. Skipping overwrite.")
        return

    with open(tsv_file, "a") as f:  # Open in append mode
        if os.stat(tsv_file).st_size == 0:  # Check if file is empty
            header = list(metrics[0].keys())
            f.write("\t".join(header) + "\n")  # Write header only if file is empty
        for metric in metrics:
            line = "\t".join([str(val) for val in metric.values()])
            f.write(line + "\n")


def process_reads_bam(bam_file, tsv_file, tags, chunk_size=10000, force_overwrite=True):
    print(f"{time.strftime('%D - %H:%M:%S', time.localtime())} | Processing reads...")
    for metrics_chunk in calculate_metrics_bam(bam_file, tags, chunk_size):
        print(f"Processed {len(metrics_chunk)} reads...")
        print(
            f"{time.strftime('%D - %H:%M:%S', time.localtime())} | Writing metrics to file..."
        )
        write_tsv(metrics_chunk, tsv_file, force_overwrite)


def main_bam(bam_file, tsv_file, tags, chunk_size=10000, force_overwrite=True):
    output_dir = os.path.dirname(tsv_file)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    process_reads_bam(bam_file, tsv_file, tags, chunk_size, force_overwrite)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Parse BAM/SAM file and write metrics to TSV."
    )
    parser.add_argument("--bam_file", type=str, help="Path to the input BAM/SAM file.")
    parser.add_argument("--tsv_file", type=str, help="Path to the output TSV file.")
    parser.add_argument(
        "--tags",
        type=str,
        nargs="+",
        default=[],
        help="List of tags to include in the output TSV. Default is an empty list.",
    )
    parser.add_argument(
        "--chunk-size",
        type=int,
        default=10000,
        help="Number of reads to process in each chunk. Default is 10000.",
    )
    parser.add_argument(
        "--force-overwrite",
        action="store_true",
        default=True,
        help="Force overwrite of the output file if it exists. Default is True.",
    )
    args = parser.parse_args()

    main_bam(
        args.bam_file, args.tsv_file, args.tags, args.chunk_size, args.force_overwrite
    )
