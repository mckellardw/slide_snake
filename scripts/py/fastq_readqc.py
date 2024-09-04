import sys
import argparse
import gzip
import os
from itertools import groupby
from concurrent.futures import ProcessPoolExecutor, as_completed
import time

# Usage:
"""
python scripts/py/fastq_readqc.py \
    data/test_seeker/heart_ctrl_4k_ont.fq.gz \
    heart_ctrl_4k_ont.tsv \
    --cores 4 \
    --chunk_size 1000
"""

def count_reads_in_fastq(fastq_file):
    """
    Count the number of reads in a gzipped FASTQ file.
    
    Args:
    fastq_file (str): Path to the gzipped FASTQ file.
    
    Returns:
    int: Number of reads in the FASTQ file.
    """
    is_compressed = fastq_file.endswith(".gz")
    line_count = 0
    
    with gzip.open(fastq_file, "rt") if is_compressed else open(fastq_file, "r") as f:
        # Can't use '@' to count lines in these b/c ONT used '@' in their Q scores...... DUMB
        for line in f:
            line_count += 1
    
    read_count = line_count // 4
    return read_count

def remove_file_if_exists(file_path):
    """
    Check if a file exists and remove it if it does.

    Args:
    file_path (str): Path to the file to be checked and potentially removed.

    Returns:
    bool: True if the file was removed, False if it didn't exist.
    """
    if os.path.exists(file_path):
        try:
            os.remove(file_path)
            print(f"{time.strftime('%D - %H:%M:%S', time.localtime())} | Removed old output file [{file_path}]...")
            return True
        except OSError as e:
            print(f"{time.strftime('%D - %H:%M:%S', time.localtime())} | Error: {e.strerror}. Unable to remove '{file_path}'.")
            return False
    else:
        return False

def calculate_metrics(fastq_file, start, chunk_size=100000):
    is_compressed = fastq_file.endswith(".gz")
    metrics = []
    read_id = "TMP"
    offset = start * 4  # Adjust for 4 lines per read
    # print(f"Offset: {offset}")

    with gzip.open(fastq_file, "rt") if is_compressed else open(fastq_file, "r") as f:
        # f.seek(offset) # start at offset
        for i, line in enumerate(f):
            # Skip to offset
            if i >= offset:
                if i % 4 == 0 and line.startswith("@"):  # ID line
                    read_id = line.strip().replace("@", "").split(" ", 1)[0]
                
                if i % 4 == 1:  # Sequence line
                    seq = line.strip()
                    if len(seq) > 0:
                        first_base = seq[0]
                        last_base = seq[-1]
                        gc_count = seq.count("G") + seq.count("C")
                        gc_percent = round(gc_count * 100 / len(seq), 2)
                        purine_count = seq.count("G") + seq.count("A")
                        purine_percent = round(purine_count * 100 / len(seq), 2)
                        homopolymer_lengths = [len(list(g)) for k, g in groupby(seq)]
                        longest_homopolymer = max(homopolymer_lengths)
                        homopolymer_base = seq[
                            homopolymer_lengths.index(longest_homopolymer)
                        ]
                    else:
                        first_base = None
                        last_base = None
                        gc_percent = None
                        purine_percent = None
                        homopolymer_lengths = None
                        longest_homopolymer = None
                        homopolymer_base = None
                    # elif i % 4 == 3:  # Quality line
                    #TODO

                    metrics.append(
                        {
                            "Read_ID": read_id,
                            "Read_Length": len(seq),
                            "GC_Percent": gc_percent,
                            "Purine_Percent": purine_percent,
                            "First_Base": first_base,
                            "Last_Base": last_base,
                            "Longest_Homopolymer": longest_homopolymer,
                            "Homopolymer_Base": homopolymer_base,
                        }
                    )

                #TODO- use all cores, not just 1 chunk per core...
                # if (i + 1) % chunk_size == 0: 
                    # yield metrics
                    # metrics = []
                
                # Exit after 1 chunk
                if i == (offset + (4*chunk_size)):
                    print(f"Offset: {offset} to {i} -> {i-offset} lines | {(i-offset)/4} reads")
                    return metrics

    if metrics:
        print(f"Last Offset: {offset} to {i} -> {i-offset} lines | {(i-offset)/4} reads")
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
    for read_metrics in calculate_metrics(fastq_file, start, chunk_size):
        write_tsv(read_metrics, tsv_file)


def process_reads(fastq_file, tsv_file, chunk_size=100000, cores=1):
    # Determine the number of chunks to process
    print(f"{time.strftime('%D - %H:%M:%S', time.localtime())} | Counting reads...")
    total_reads = count_reads_in_fastq(fastq_file)
    print(f"{time.strftime('%D - %H:%M:%S', time.localtime())} | {total_reads} total reads found...")

    total_chunks = (
        total_reads + chunk_size - 1
    ) // chunk_size  # Round up to cover all reads
    print(f"{time.strftime('%D - %H:%M:%S', time.localtime())} | Calculating across {total_chunks} chunks...")

    with ProcessPoolExecutor(max_workers=cores) as executor:
        futures = [
            executor.submit(worker, fastq_file, tsv_file, i * chunk_size, chunk_size)
            for i in range(total_chunks)
        ]

        completed_chunks = 0
        for future in as_completed(futures):
            try:
                result = future.result()
                completed_chunks += 1
                print(f"{time.strftime('%D - %H:%M:%S', time.localtime())} | Completed chunks {completed_chunks}/{total_chunks}")
            except Exception as e:
                print(f"{time.strftime('%D - %H:%M:%S', time.localtime())} | Error in chunk: {str(e)}")



def main(fastq_file, tsv_file, cores, chunk_size):
    # Make output directory if needed
    output_dir = os.path.dirname(tsv_file)
    if len(output_dir) > 0 and not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    output_removed = remove_file_if_exists(tsv_file)
        
    print(f"{time.strftime('%D - %H:%M:%S', time.localtime())} | Processing reads across {cores} cores...")

    process_reads(fastq_file, tsv_file, chunk_size, cores)
    
    print(f"{time.strftime('%D - %H:%M:%S', time.localtime())} | Done!")


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

    print(
        f"Input fastq:      {args.fastq_file}\n"
        f"Output tsv:       {args.tsv_file}\n"
        f"Number of cores:  {args.cores}\n"
        f"Chunk size:       {args.chunk_size}\n"
        f"\n"
    )

    # if args.cores == 1:

    # elif args.cores > 1:
    main(args.fastq_file, args.tsv_file, args.cores, args.chunk_size)
