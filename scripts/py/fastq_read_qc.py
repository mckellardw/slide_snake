import sys
import argparse
import gzip
import os
from itertools import groupby
from concurrent.futures import ProcessPoolExecutor, as_completed
import time

# Usage:
"""
python scripts/py/fastq_read_qc2.py \
    data/test_seeker/heart_ctrl_4k_ont.fq.gz \
    heart_ctrl_4k_ont.tsv \
    --cores 4 \
    --chunk_size 1000
"""


def calculate_metrics(fastq_file, start, chunk_size=100000):
    is_compressed = fastq_file.endswith(".gz")
    metrics = []
    read_id = ""
    offset = start * 4  # Adjust for 4 lines per read

    with gzip.open(fastq_file, "rt") if is_compressed else open(fastq_file, "r") as f:
        f.seek(offset)
        for i, line in enumerate(f):
            if i % 4 == 0:  # ID line
                read_id = line.strip().replace("@", "").split(" ", 1)[0]
            elif i % 4 == 1:  # Sequence line
                seq = line.strip()
                if len(seq) > 0:
                    first_base = seq[0]
                    last_base = seq[-1]
                    gc_count = seq.count("G") + seq.count("C")
                    gc_percent = round(gc_count * 100 / len(seq), 2)
                    homopolymer_lengths = [len(list(g)) for k, g in groupby(seq)]
                    longest_homopolymer = max(homopolymer_lengths)
                    homopolymer_base = seq[
                        homopolymer_lengths.index(longest_homopolymer)
                    ]
                else:
                    first_base = None
                    last_base = None
                    gc_percent = None
                    homopolymer_lengths = None
                    longest_homopolymer = None
                    homopolymer_base = None

                metrics.append(
                    {
                        "Read_ID": read_id,
                        "Read_Length": len(seq),
                        "GC_Percent": gc_percent,
                        "First_Base": first_base,
                        "Last_Base": last_base,
                        "Longest_Homopolymer": longest_homopolymer,
                        "Homopolymer_Base": homopolymer_base,
                    }
                )

            if (i + 1) % chunk_size == 0:
                yield metrics
                metrics = []

        if metrics:
            yield metrics


def write_tsv(metrics, tsv_file):
    file_is_empty = not os.path.exists(tsv_file) or os.stat(tsv_file).st_size == 0

    with open(tsv_file, "a") as f:
        if file_is_empty:
            f.write("\t".join(metrics[0].keys()) + "\n")
        for metric in metrics:
            line = "\t".join([str(val) for val in metric.values()])
            f.write(line + "\n")


def worker(fastq_file, tsv_file, start, chunk_size):
    for metrics_chunk in calculate_metrics(fastq_file, start, chunk_size):
        print(f"Processed {len(metrics_chunk)} reads...")
        print(
            f"{time.strftime('%D - %H:%M:%S', time.localtime())} | Writing metrics to file..."
        )
        write_tsv(metrics_chunk, tsv_file)


def process_reads(fastq_file, tsv_file, chunk_size=100000, cores=1):
    print(f"{time.strftime('%D - %H:%M:%S', time.localtime())} | Processing reads...")

    # Determine the number of chunks to process
    total_reads = sum(
        1 for _ in calculate_metrics(fastq_file, 0, chunk_size=1)
    )  # Only counting reads
    total_chunks = (
        total_reads + chunk_size - 1
    ) // chunk_size  # Round up to cover all reads

    with ProcessPoolExecutor(max_workers=cores) as executor:
        futures = [
            executor.submit(worker, fastq_file, tsv_file, i * chunk_size, chunk_size)
            for i in range(total_chunks)
        ]

        for future in as_completed(futures):
            future.result()


def main(fastq_file, tsv_file, cores, chunk_size):
    output_dir = os.path.dirname(tsv_file)
    if len(output_dir) > 0 and not os.path.exists(output_dir):
        os.makedirs(output_dir)

    process_reads(fastq_file, tsv_file, chunk_size, cores)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Parse FASTQ file and write metrics to TSV."
    )
    parser.add_argument("fastq_file", type=str, help="Path to the input FASTQ file.")
    parser.add_argument("tsv_file", type=str, help="Path to the output TSV file.")
    parser.add_argument(
        "--cores",
        type=int,
        default=1,
        help="Number of cores to use for multithreading. Default is 1.",
    )
    parser.add_argument(
        "--chunk_size",
        type=int,
        default=100000,
        help="Number of reads to process in each chunk. Default is 100k.",
    )
    args = parser.parse_args()

    main(args.fastq_file, args.tsv_file, args.cores, args.chunk_size)
