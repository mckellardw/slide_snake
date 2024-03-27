import h5py
import sys
import numpy as np
from concurrent.futures import ThreadPoolExecutor


def seqDecode(seqInt, seqLen):
    # ATCG_BASES = ['A', 'C', 'T', 'G']
    ATCG_BASES = ["T", "G", "A", "C"]  # use complement to decode
    seqs = ""
    for i in range(seqLen):
        tint = (seqInt >> (i * 2)) & 3
        seqs += ATCG_BASES[tint]
    return seqs


def process_chunk(chunk, progress, seq_length, verbose):
    lines_per_update = 20000000
    results = []
    for row_idx, row in enumerate(chunk):
        for col_idx, item in enumerate(row):
            int_item = int(item)  # Convert item to integer if necessary
            if int_item == 0:  # Skip if the integer is zero
                continue
            nucleotide_pairs = seqDecode(int_item, seq_length)
            results.append(f"{nucleotide_pairs}\t{row_idx}\t{col_idx}\n")
            progress[0] += 1
            if verbose and progress[0] % lines_per_update == 0:
                print(f"Processed {progress[0]} lines")
    return results


def process_group(group, output_file, seq_length, n_threads, verbose):
    for key, value in group.items():
        if isinstance(value, h5py.Dataset):
            dataset_array = np.array(value)  # Convert the dataset to a NumPy array
            chunk_size = dataset_array.shape[0] // n_threads
            chunks = [
                dataset_array[i : i + chunk_size]
                for i in range(0, dataset_array.shape[0], chunk_size)
            ]
            progress = [0]
            results = []

            with ThreadPoolExecutor(max_workers=n_threads) as executor:
                futures = [
                    executor.submit(process_chunk, chunk, progress, seq_length, verbose)
                    for chunk in chunks
                ]
                for future in as_completed(futures):
                    results.extend(future.result())

            with open(output_file, "a") as out:
                for result in results:
                    out.write(result)

        elif isinstance(value, h5py.Group):
            process_group(value, output_file, seq_length, n_threads, verbose)


def main(input_file, output_file, seq_length=25, n_threads=16, verbose=False):
    with h5py.File(input_file, "r") as f:
        process_group(f, output_file, seq_length, n_threads, verbose)


if __name__ == "__main__":
    if len(sys.argv) < 3:
        print(
            "Usage: python h5_to_list.py input_file output_file [seq_length] [n_threads] [--verbose]"
        )
        sys.exit(1)

    input_file = sys.argv[1]
    output_file = sys.argv[2]
    seq_length = int(sys.argv[3]) if len(sys.argv) > 3 else 25
    n_threads = int(sys.argv[4]) if len(sys.argv) > 4 else 16
    verbose = "--verbose" in sys.argv

    main(input_file, output_file, seq_length, n_threads, verbose)
