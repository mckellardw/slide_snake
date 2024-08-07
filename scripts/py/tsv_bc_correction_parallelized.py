import sys
import os
import random
import string
import csv
import argparse
from itertools import product
import editdistance as ed
import sys
import csv
import argparse

# Usage:
# SlideSeq
""" 
python scripts/py/tsv_bc_correction_parallelized.py \
    --tsv_in barcodes.tsv \
    --tsv_out_full barcodes_corrected_full.tsv \
    --tsv_out_slim barcodes_corrected.tsv \
    --id_column 0 \
    --bc_columns 1 2 \
    --concat_bcs True \
    --whitelist_files out/Heart_Control/bc/whitelist.txt \
    --max_ham 2 \
    --threads 10
"""

## Visium
"""
python scripts/py/tsv_bc_correction.py \
    --tsv_in barcodes.tsv \
    --tsv_out_full barcodes_corrected_full.tsv \
    --tsv_out_slim barcodes_corrected.tsv \
    --id_column 0 \
    --bc_columns 1 \
    --concat_bcs True \
    --whitelist_files out/Vis_yPAP_3D/bc/whitelist.txt \
    --max_ham 2 \
    --threads 10
"""


# Random name for temp directory
def generate_temp_dir_name(length=16):
    # Combine all ASCII letters and digits
    chars = string.ascii_letters + string.digits
    # Generate a random string of the specified length
    return "temp_" + "".join(random.choice(chars) for _ in range(length))


# Get number of lines in a file
def file_line_count(filename):
    with open(filename, "r") as f:
        for i, _ in enumerate(f):
            pass
    return i + 1


# Split a txt file into chunks for parallelization
def split_file(file_path, temp_dir, num_chunks):
    # Ensure the temporary directory exists
    if not os.path.exists(temp_dir):
        os.makedirs(temp_dir)

    # Calculate the size of each chunk (number of lines)
    chunk_size = (file_line_count(file_path) // num_chunks) + 1
    # print(chunk_size)

    # Initialize an empty list to store chunk file names
    chunk_file_names = [
        f"{temp_dir}/chunk_{str(i).zfill(3)}.tsv" for i in range(num_chunks)
    ]

    # Open the original file and start writing the chunks
    with open(file_path, "r") as infile:
        chunk_index = 0
        line_count = 0
        for line in infile:
            # Determine the output file name
            output_file_name = chunk_file_names[chunk_index]

            # Open the output file
            with open(output_file_name, "a") as outfile:
                outfile.write(line)

            line_count += 1

            # Check if we have reached the next chunk
            if line_count >= chunk_size:
                line_count = 0
                chunk_index += 1

    return chunk_file_names


# Concatenate a list of files
def concatenate_files(filenames, output_filename):
    """
    Concatenates a list of files into a single output file.

    Parameters:
    - filenames: A list of file names to concatenate.
    - output_filename: The name of the output file.
    """
    # Open the output file in write mode
    with open(output_filename, "w") as outfile:
        # Iterate through the list of input files
        for filename in filenames:
            # Open each input file in read mode
            with open(filename, "r") as infile:
                # Read the content of the input file
                content = infile.read()
                # Write the content to the output file
                outfile.write(content)


# Simple python version of R rep() function
def rep(val, n):
    return [val for i in range(0, n)]


def filter_whitelist_by_kmers(wl, kmers, kmer_to_bc_index):
    """
    Given a list of whitelisted barcodes, return just the
    subset that contain any of the kmers contained in the
    query barcode.

    :param wl: Full barcode whitelist
    :type wl: list
    :param kmers: K-mers to use for whitelist filtering
    :type kmers: list
    :param kmer_to_bc_index: Map of k-mers to the whitelist indices corresponding
        to all barcodes containing that k-mer
    :type kmer_to_bc_index: dict
    :return: List of filtered barcodes
    :rtype: list
    """
    # collect sets of indices that each kmer points to
    id_sets = [
        kmer_to_bc_index[kmer] for kmer in kmers if kmer in kmer_to_bc_index.keys()
    ]

    # retain all barcodes that have at least one kmer match with the query barcode
    all_filt_indices = list(set().union(*id_sets))
    filt_wl = [wl[i] for i in all_filt_indices]
    return filt_wl


def split_seq_into_kmers(seq, k):
    """
    Decompose the supplied <seq> into N=len(seq)-k+1 k-mers.

    :param seq: String of nucleotides
    :type seq: str
    :param k: k-mer length
    :type k: int
    :return: List of k-mers
    :rtype: list
    """
    assert len(seq) >= k, "Pick a value for k that is less than len(barcode)"

    kmers = []
    for i in range(0, len(seq) - k + 1):
        kmer = seq[i : i + k]
        kmers.append(kmer)
    return kmers


def load_whitelist(whitelist, k=5):
    """
    Read in barcode whitelist and create dictionary mapping each k-mer to all
    barcodes in the whitelist containing that k-mer.

    :param whitelist: Path to the barcode whitelist
    :type whitelist: str
    :param k: k-mer length
    :type k: int
    :return: List of whitelisted barcodes and dictionary mapping all k-mers to
        indices in the whitelist corresponding to barcodes containing that k-mer
    :rtype: list, dict
    """
    wl = []
    with open(whitelist) as file:
        for line in file:
            bc = line.strip().split("-")[0]
            wl.append(bc)

    wl.sort()
    kmer_to_bc_index = {}
    for index, bc in enumerate(wl):
        bc_kmers = split_seq_into_kmers(bc, k)
        for bc_kmer in bc_kmers:
            if bc_kmer not in kmer_to_bc_index.keys():
                kmer_to_bc_index[bc_kmer] = set([index])
            else:
                kmer_to_bc_index[bc_kmer].add(index)
    return wl, kmer_to_bc_index


def calc_ed_with_whitelist(bc_uncorr, whitelist):
    """
    Find minimum and runner-up barcode edit distance by iterating through the
    whitelist of expected barcodes.

    :param bc_uncorr: Uncorrected cell barcode
    :type bc_uncorr: str
    :param whitelist: Filtered whitelist of cell barcodes
    :type whitelist: list
    :return: Corrected barcode assignment, edit distance, and difference in edit
        distance between the top match and the next closest match
    :rtype: str, int, int
    """
    bc_match = "X" * len(bc_uncorr)
    bc_match_ed = len(bc_uncorr)
    next_bc_match_ed = len(bc_uncorr)
    for wl_bc in whitelist:
        d = ed.eval(bc_uncorr, wl_bc)  # Use the ed module here
        if d < bc_match_ed:
            next_bc_match_ed = bc_match_ed
            bc_match_ed = d
            bc_match = wl_bc
        elif d < next_bc_match_ed:
            next_bc_match_ed = d
    next_match_diff = next_bc_match_ed - bc_match_ed

    return bc_match, bc_match_ed, next_match_diff


def process_tsv(
    tsv_in,
    tsv_out_full,
    tsv_out_slim,
    id_column,
    bc_columns,
    concat_bcs,
    whitelist_files,
    max_hams,
    verbose=False,
):
    # prep whitelists
    whitelists = {}
    for i, whitelist_file in enumerate(whitelist_files):
        wl, kmer_to_bc_index = load_whitelist(
            whitelist_file, k=5
        )  # Assuming k=5 for simplicity
        whitelists[i] = wl

    # Prepare the output file
    with open(tsv_out_full, "w", newline="") as outfile_full:
        writer_full = csv.writer(outfile_full, delimiter="\t")
        # writer.writerow(['Read_ID', 'Original_Barcode', 'Corrected_Barcode', 'Hamming_Distance'])
        with open(tsv_out_slim, "w", newline="") as outfile_slim:
            writer_slim = csv.writer(outfile_slim, delimiter="\t")
            # writer.writerow(['Read_ID    Corrected_Barcode'])

            read_count = 0
            with open(tsv_in, "r") as infile:
                reader = csv.reader(infile, delimiter="\t")
                next(reader)  # Skip header row

                for row in reader:
                    read_count += 1
                    if verbose and read_count % 1000000 == 0:
                        print(f"{read_count} reads processed...")

                    read_id = row[id_column]
                    barcodes = [row[i] for i in bc_columns]

                    if concat_bcs:
                        barcodes = ["".join(barcodes)]

                    row2write = [read_id]  # initialize row to write in output .tsv

                    RANGE = [
                        [barcodes[i], whitelists[i], max_hams[i]]
                        for i in range(len(barcodes))
                    ]
                    for barcode, whitelist, max_ham in RANGE:
                        (
                            corrected_bc,
                            bc_match_ham,
                            next_match_diff,
                        ) = calc_ed_with_whitelist(barcode, whitelist)

                        if bc_match_ham <= max_ham:
                            # row2write.append(f"{barcode}\t{corrected_bc}\t{bc_match_ed}")
                            row2write.extend(
                                [barcode, corrected_bc, bc_match_ham, next_match_diff]
                            )
                        else:
                            # row2write.append(f"{barcode}\t \t{bc_match_ed}")
                            row2write.extend(
                                [barcode, "-", bc_match_ham, next_match_diff]
                            )
                    writer_full.writerow(row2write)
                    writer_slim.writerow([read_id, corrected_bc])


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Correct barcodes in a TSV file.")
    parser.add_argument("--tsv_in", required=True, help="Path to the input TSV file.")
    parser.add_argument(
        "--tsv_out_full",
        required=True,
        help="Path to the output TSV file. Columns contain ['Read_ID', 'Original_Barcode', 'Corrected_Barcode', 'Hamming_Distance']",
    )
    parser.add_argument(
        "--tsv_out_slim",
        required=True,
        help="Path to output TSV file which only contains the read ID and final BC. Columns contain ['Read_ID', 'Corrected_Barcode']",
    )
    parser.add_argument(
        "--whitelist_files",
        nargs="+",
        required=True,
        help="Space-separated list of whitelist file paths.",
    )
    parser.add_argument(
        "--id_column",
        nargs="+",
        type=int,
        default=0,
        help="Column in .tsv corresponding to the read IDs (default: 0).",
    )
    parser.add_argument(
        "--bc_columns",
        nargs="+",
        type=int,
        required=True,
        help="Columns in .tsv corresponding to the uncorrected barcodes.",
    )
    parser.add_argument(
        "--concat_bcs",
        type=bool,
        default=True,
        help="Columns in .tsv corresponding to the uncorrected barcodes.",
    )
    parser.add_argument(
        "--max_hams",
        nargs="+",
        type=int,
        default=2,
        help="Maximum Hamming distance for correction (default: 2).",
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=1,
        help="Number of threads to use (default: 1).",
    )

    args = parser.parse_args()

    # param checks ------
    if len(args.max_hams) != len(args.bc_columns) and len(args.max_hams) == 1:
        args.max_hams = rep(val=args.max_hams[0], n=len(args.bc_columns))

    if args.concat_bcs and len(args.whitelist_files) > 1:
        print(f"Need a merged barcode whitelist!")
        sys.exit(1)
    
    # Print run settings for log files ----
    print(
        f"Input tsv:                    {args.tsv_in}\n"
        f"Output tsv (Full Info):       {args.tsv_out_full}\n"
        f"Output tsv (CBs only):        {args.tsv_out_slim}\n"
        f"Read ID column:               {args.id_column}\n"
        f"Barcode column(s):            {args.bc_columns}\n"
        f"Concatenate uncorrected BCs?: {args.concat_bcs}\n"
        f"Whitelist file(s):            {args.whitelist_files}\n"
        f"Maximum hamming distance(s):  {args.max_hams}\n"
    )

    # Single-threaded = verbose
    if args.threads == 1:
        print(f"Running on {args.threads} thread...")
        print("")
        process_tsv(
            tsv_in=args.tsv_in,
            tsv_out_full=args.tsv_out_full,
            tsv_out_slim=args.tsv_out_slim,
            # id_column=args.id_column[0],
            id_column=args.id_column,
            bc_columns=args.bc_columns,
            concat_bcs=args.concat_bcs,
            whitelist_files=args.whitelist_files,
            max_hams=args.max_hams,
            verbose=True,
        )

    # Multi-threaded = not verbose
    elif args.threads > 1:
        print(f"Running on {args.threads} threads...")
        # Source: https://superfastpython.com/multiprocessing-pool-for-loop/
        import multiprocessing

        print("")
        print(f"Splitting input tsv...")
        temp_dir = generate_temp_dir_name()
        print(f"Temporary directory: {temp_dir}")
        temp_tsvs_in = split_file(
            file_path=args.tsv_in, temp_dir=temp_dir, num_chunks=args.threads
        )

        print(f"Correcting barcodes...")

        temp_tsvs_out_full = [
            fn.replace(".tsv", "_corr_full.tsv") for fn in temp_tsvs_in
        ]
        temp_tsvs_out_slim = [
            fn.replace(".tsv", "_corr_slim.tsv") for fn in temp_tsvs_in
        ]

        # multiprocessing
        items = [
            (
                temp_tsvs_in[i],
                temp_tsvs_out_full[i],
                temp_tsvs_out_slim[i],
                # args.id_column[0],
                args.id_column,
                args.bc_columns,
                args.concat_bcs,
                args.whitelist_files,
                args.max_hams,
            )
            for i in list(range(args.threads))
        ]

        with multiprocessing.Pool(args.threads) as pool:
            multi_out = pool.starmap(process_tsv, items)

        # process_tsv(
        #     tsv_in=args.tsv_in,
        #     tsv_out_full=args.tsv_out_full,
        #     tsv_out_slim=args.tsv_out_slim,
        #     id_column=args.id_column[0],
        #     bc_columns=args.bc_columns,
        #     concat_bcs=args.concat_bcs,
        #     whitelist_files=args.whitelist_files,
        #     max_hams=args.max_hams
        # )

        print(f"Concatenating results...")
        concatenate_files(
            filenames=temp_tsvs_out_full, output_filename=args.tsv_out_full
        )
        concatenate_files(
            filenames=temp_tsvs_out_slim, output_filename=args.tsv_out_slim
        )

        print(f"Removing temp files (temp_dir: `{temp_dir}`)...")
        os.system(f"rm -rf {temp_dir}") if os.path.exists(temp_dir) else None

        print(f"Done!")
    else:
        print(f"Incorrect number of threads specified [{args.threads}]...")
