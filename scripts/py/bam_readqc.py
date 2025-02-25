import argparse
import os
import pysam
import time
from itertools import groupby
from collections import Counter
from concurrent.futures import ProcessPoolExecutor, as_completed

# Example usage:
"""
python scripts/py/bam_read_qc.py \
    --tags {params.TAGS} \
    --chunk-size {params.CHUNK_SIZE} \
    --bam_file {input.BAM} \
    --tsv_file {output.TSV} 
"""


def calculate_metrics_bam(bam_file, tags, chunk_start, chunk_size=10000):
    metrics = []
    with pysam.AlignmentFile(bam_file, "rb") as samfile:
        for i, read in enumerate(samfile):
            if i < chunk_start:
                continue
            if i >= chunk_start + chunk_size:
                break
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


def worker_bam(bam_file, tags, chunk_start, chunk_size, temp_dir):
    temp_file = os.path.join(temp_dir, f"chunk_{chunk_start}.tsv")
    try:
        for metrics_chunk in calculate_metrics_bam(bam_file, tags, chunk_start, chunk_size):
            write_tsv(metrics_chunk, temp_file)
        print(
            f"{time.strftime('%D - %H:%M:%S', time.localtime())} | Worker completed for chunk starting at {chunk_start}"
        )
    except Exception as e:
        print(
            f"{time.strftime('%D - %H:%M:%S', time.localtime())} | Error in worker for chunk starting at {chunk_start}: {str(e)}"
        )
        raise

def merge_tsv_files(temp_dir, final_tsv_file):
    if not os.path.exists(temp_dir):
        raise FileNotFoundError(f"Temporary directory {temp_dir} does not exist.")
    with open(final_tsv_file, "w") as outfile:
        for i, temp_file in enumerate(sorted(os.listdir(temp_dir))):
            temp_file_path = os.path.join(temp_dir, temp_file)
            with open(temp_file_path, "r") as infile:
                if i == 0:
                    outfile.write(infile.read())
                else:
                    next(infile)  # Skip header
                    outfile.write(infile.read())

def process_reads_bam(bam_file, tsv_file, tags, chunk_size=10000, cores=1, force_overwrite=True):
    print(f"{time.strftime('%D - %H:%M:%S', time.localtime())} | Processing reads...")
    temp_dir = os.path.join(os.path.dirname(tsv_file), "temp_chunks")
    os.makedirs(temp_dir, exist_ok=True)

    with pysam.AlignmentFile(bam_file, "rb") as samfile:
        total_reads = sum(1 for _ in samfile)
    total_chunks = (total_reads + chunk_size - 1) // chunk_size

    with ProcessPoolExecutor(max_workers=cores) as executor:
        futures = [
            executor.submit(worker_bam, bam_file, tags, i * chunk_size, chunk_size, temp_dir)
            for i in range(total_chunks)
        ]

        completed_chunks = 0
        for future in as_completed(futures):
            try:
                result = future.result()
                completed_chunks += 1
                print(
                    f"{time.strftime('%D - %H:%M:%S', time.localtime())} | Completed chunks {completed_chunks}/{total_chunks} ({(completed_chunks/total_chunks)*100:.2f}%)"
                )
            except Exception as e:
                print(
                    f"{time.strftime('%D - %H:%M:%S', time.localtime())} | Error in chunk: {str(e)}"
                )

    merge_tsv_files(temp_dir, tsv_file)
    print(
        f"{time.strftime('%D - %H:%M:%S', time.localtime())} | Merged all chunks into {tsv_file}"
    )
    print(f"{time.strftime('%D - %H:%M:%S', time.localtime())} | Total reads processed: {total_reads}")
    # Clean up temporary directory
    # for temp_file in os.listdir(temp_dir):
    #     os.remove(os.path.join(temp_dir, temp_file))
    # os.rmdir(temp_dir)

def parse_args():
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
        "--cores",
        type=int,
        default=1,
        help="Number of cores to use for multithreading. Default is 1.",
    )
    parser.add_argument(
        "--force-overwrite",
        action="store_true",
        default=True,
        help="Force overwrite of the output file if it exists. Default is True.",
    )
    return parser.parse_args()

def main_bam(bam_file, tsv_file, tags, chunk_size=10000, cores=1, force_overwrite=True):
    output_dir = os.path.dirname(tsv_file)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    process_reads_bam(bam_file, tsv_file, tags, chunk_size, cores, force_overwrite)

if __name__ == "__main__":
    args = parse_args()

    print(
        f"Input BAM:        {args.bam_file}\n"
        f"Output TSV:       {args.tsv_file}\n"
        f"Tags:             {args.tags}\n"
        f"Chunk size:       {args.chunk_size}\n"
        f"Number of cores:  {args.cores}\n"
        f"Force overwrite:  {args.force_overwrite}\n"
        f"\n"
    )

    main_bam(
        args.bam_file, args.tsv_file, args.tags, args.chunk_size, args.cores, args.force_overwrite
    )
