import sys
import argparse
import gzip
import os
from collections import Counter
from itertools import groupby
from concurrent.futures import ThreadPoolExecutor, as_completed

import time


def calculate_metrics(fastq_file, individual=True):
    metrics = []
    read_ids = []
    read_id = ""

    is_compressed = fastq_file.endswith(".gz")

    with gzip.open(fastq_file, "rt") if is_compressed else open(fastq_file, "r") as f:
        for i, line in enumerate(f):
            if i % 4 == 0:  # ID line
                # read_ids.append(line.strip())
                read_id = line.strip()
            elif i % 4 == 1:  # Sequence line
                seq = line.strip()
                if len(seq) > 0:
                    first_base = seq[0]
                    last_base = seq[-1]
                    gc_count = seq.count("G") + seq.count("C")
                    gc_percent = round(gc_count * 100 / len(seq), 2)
                    homopolymer_lengths = [len(list(g)) for k, g in groupby(seq)]
                    longest_homopolymer = max(homopolymer_lengths)
                    homopolymer_base = seq[homopolymer_lengths.index(longest_homopolymer)]
                else:                    
                    first_base = None
                    last_base = None
                    gc_percent = None
                    homopolymer_lengths = None
                    longest_homopolymer = None
                    homopolymer_base = None
                
                
                metrics.append(
                    {
                        # "Read_ID": read_ids[-1],
                        "Read_ID": read_id,
                        "Read_Length": len(seq),
                        "GC_Percent": gc_percent,
                        "First_Base": first_base,
                        "Last_Base": last_base,
                        "Longest_Homopolymer": longest_homopolymer,
                        "Homopolymer_Base": homopolymer_base,
                    }
                )

    return metrics


def write_tsv(metrics, tsv_file):
    with open(tsv_file, "w") as f:
        f.write(
            "\t".join(metrics[0].keys())+"\n"
        )
        for metric in metrics:
            line="\t".join([str(val) for val in metric.values()])
            f.write(line+"\n")


def process_reads(fastq_file, tsv_file):
    print(f"{time.strftime('%D - %H:%M:%S', time.localtime())} | Processing reads...")
    metrics = calculate_metrics(fastq_file)
    print(f"Processed {len(metrics)} reads...")
    print(f"{time.strftime('%D - %H:%M:%S', time.localtime())} | Writing metrics to file...")
    write_tsv(metrics, tsv_file)


# def main(fastq_file, tsv_file, threads):
#     output_dir = os.path.dirname(tsv_file)
#     if not os.path.exists(output_dir):
#         os.makedirs(output_dir)

#     with ThreadPoolExecutor(max_workers=threads) as executor:
#         executor.submit(process_reads, fastq_file, tsv_file)

def main(fastq_file, tsv_file, threads):
    output_dir = os.path.dirname(tsv_file)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # if threads > 1:
    #     with ThreadPoolExecutor(max_workers=threads) as executor:
    #         futures = [executor.submit(process_reads, fastq_file, tsv_file)]
    #         for future in as_completed(futures):
    #             pass
    # else:
    #     process_reads(fastq_file, tsv_file)
    process_reads(fastq_file, tsv_file)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Parse FASTQ file and write metrics to TSV."
    )
    parser.add_argument(
        "fastq_file", 
        type=str, 
        help="Path to the input FASTQ file."
    )
    parser.add_argument(
        "tsv_file", 
        type=str, 
        help="Path to the output TSV file."
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=1,
        help="Number of threads to use for multithreading. Default is 1.",
    )
    args = parser.parse_args()

    main(args.fastq_file, args.tsv_file, args.threads)
