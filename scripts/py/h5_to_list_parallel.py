import h5py
import sys
import numpy as np
from concurrent.futures import ThreadPoolExecutor

def seqDecode(seqInt, seqLen):
    ATCG_BASES = ['A', 'C', 'T', 'G']
    seqs = ""
    for i in range(seqLen):
        tint = (seqInt >> (i * 2)) & 3
        seqs += ATCG_BASES[tint]
    return seqs

def process_chunk(chunk, progress, seq_length, verbose):
    lines_per_update = 10000000
    results = []
    for row_idx, row in enumerate(chunk):
        for col_idx, item in enumerate(row):
            int_item = int(item)  # Convert item to integer if necessary
            nucleotide_pairs = seqDecode(int_item, seq_length)
            results.append(f"{nucleotide_pairs}\t{row_idx}\t{col_idx}\n")
            progress[0] += 1
            if verbose and progress[0] % lines_per_update == 0:
                print(f"Processed {progress[0]} lines")
    return results


def process_group(group, output_file, seq_length=25, n_chunks=4, verbose=False):
    progress = [0]  # Using a list to allow the integer to be mutable inside the process_chunk function

    for key, value in group.items():
        if isinstance(value, h5py.Dataset):
            dataset_array = np.array(value)
            chunk_size = len(dataset_array) // n_chunks

            with ThreadPoolExecutor() as executor:
                futures = []
                for chunk_start in range(0, len(dataset_array), chunk_size):
                    chunk_end = min(chunk_start + chunk_size, len(dataset_array))
                    chunk = dataset_array[chunk_start:chunk_end]
                    futures.append(executor.submit(process_chunk, chunk, progress, seq_length, verbose))

                for future in futures:
                    chunk_results = future.result()
                    for result in chunk_results:
                        output_file.write(result)

        elif isinstance(value, h5py.Group):
            process_group(value, output_file, seq_length, n_chunks, verbose)

def main(input_file, output_file, seq_length=25, n_chunks=4, verbose=False):
    with h5py.File(input_file, "r") as f, open(output_file, "w") as out:
        process_group(f, out, seq_length, n_chunks, verbose)

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python h5_to_list.py input_file output_file [--seq-length N] [--chunks N] [--verbose]")
        sys.exit(1)

    seq_length = 25
    n_chunks = 4
    verbose = False

    for idx, arg in enumerate(sys.argv):
        if arg == '--seq-length' and idx + 1 < len(sys.argv):
            seq_length = int(sys.argv[idx + 1])
        elif arg == '--chunks' and idx + 1 < len(sys.argv):
            n_chunks = int(sys.argv[idx + 1])
        elif arg == '--verbose':
            verbose = True

    main(sys.argv[1], sys.argv[2], seq_length, n_chunks, verbose)
