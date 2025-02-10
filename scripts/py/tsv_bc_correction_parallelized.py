import sys
import os
import csv
import random
import string
import argparse

# from itertools import product
# import editdistance as ed
from Levenshtein import distance, hamming
import time

startTime = time.time()

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
    --max_levens 3 \
    --min_next_match_diffs 1 \
    --k 7 \
    --threads 56
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
    --max_levens 3 \
    --min_next_match_diffs 1 \
    --k 7 \
    --threads 56
"""

## microST
"""
python scripts/py/tsv_bc_correction_parallelized.py \
    --tsv_in out/E109_A1_pool_txg/ont/barcodes_umis/microST_ssv1/read_barcodes_filtered.tsv \
    --tsv_out_full out/E109_A1_pool_txg/ont/barcodes_umis/microST_ssv1/read_barcodes_corrected_full.tsv \
    --tsv_out_slim out/E109_A1_pool_txg/ont/barcodes_umis/microST_ssv1/read_barcodes_corrected.tsv \
    --id_column 0 \
    --bc_columns 1 2 \
    --whitelist_files out/E109_A1_pool_txg/bc/whitelist_1.txt out/E109_A1_pool_txg/bc/whitelist_2.txt \
    --max_levens 3 \
    --min_next_match_diffs 1 \
    --k 5 \
    --threads 56 \
    1> out/E109_A1_pool_txg/ont/misc_logs/microST_ssv1/1c_tsv_bc_correction.log \
    2> out/E109_A1_pool_txg/ont/misc_logs/microST_ssv1/1c_tsv_bc_correction.err
"""

# fast levenshtein implementation- https://github.com/rapidfuzz/Levenshtein


def currentTime():
    return time.strftime("%D | %H:%M:%S", time.localtime())


def parse_args():
    parser = argparse.ArgumentParser(description="Correct barcodes in a TSV file.")
    parser.add_argument("--tsv_in", required=True, help="Path to the input TSV file.")
    parser.add_argument(
        "--tsv_out_full",
        required=True,
        help="Path to the output TSV file. Columns don't have titles, but contain ['Read_ID', 'Original_Barcode_1', 'Corrected_Barcode_1', 'Edit_Distance_1',, 'Original_Barcode_2', 'Corrected_Barcode_2', 'Edit_Distance_2', etc.]",
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
        action="store_true",
        help="Whether or not to combine sub-barcodes prior to correction (True for SlideSeq; False for microST/dBIT)",
    )
    parser.add_argument(
        "--max_levens",
        nargs="+",
        type=int,
        default=3,
        help="Maximum Levenshtein distance for correction; higher value will give looser correction. Space-delimited list for multiple barcode values (default: 3).",
    )
    parser.add_argument(
        "--min_next_match_diffs",
        nargs="+",
        type=int,
        default=1,
        help="Maximum difference between first and second closest matches for correction; higher value will retain more ambiguous matches. Space-delimited list for multiple barcode values (default: 1).",
    )
    parser.add_argument(
        "--k",
        type=int,
        default=5,
        help="K value to index/filter whitelist (default: 5).",
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=1,
        help="Number of threads to use (default: 1).",
    )

    args = parser.parse_args()

    return args


# Random name for temp directory
def generate_temp_dir_name(length=16):
    """
    Make a temperorary directory. Return path.
    """
    # Combine all ASCII letters and digits
    chars = string.ascii_letters + string.digits
    # Generate a random string of the specified length
    return "tmp/" + "".join(random.choice(chars) for _ in range(length))


def file_line_count(filename):
    """
    Get number of lines in a file

    Parameters:
    - filename:
    """
    with open(filename, "r") as f:
        for i, _ in enumerate(f):
            pass
    return i + 1


def split_file(file_path, temp_dir, num_chunks):
    """
    Split a txt file into chunks for parallelization.

    Parameters:
    - file_path:
    - temp_dir:
    - num_chunks:
    """
    # Ensure the temporary directory exists
    if not os.path.exists(temp_dir):
        os.makedirs(temp_dir)

    # Calculate the size of each chunk (number of lines)
    n_bcs = file_line_count(file_path)
    chunk_size = (n_bcs // num_chunks) + 1
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

    return chunk_file_names, n_bcs


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


def rep(val, n):
    """
    Simple python version of R rep() function

    Parameters:
    - val: A value to repeat.
    - n: Number of times to repeat it.
    """
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


def find_kmers(seq, k):
    """
    Decompose the supplied <seq> into N=len(seq)-k+1 k-mers.

    :param seq: String of nucleotides
    :type seq: str
    :param k: k-mer length
    :type k: int
    :return: List of k-mers
    :rtype: list
    """
    assert (
        len(seq) >= k
    ), f"Pick a value for k that is less than len(barcode); [{k}] is bigger than [{seq}]!"

    return [seq[i : i + k] for i in range(len(seq) - k + 1)]


def generate_kmer_index(whitelist, k):
    """
    Create dictionary mapping each k-mer to all
    barcodes in the whitelist containing that k-mer.

    :param whitelist: List of whitelisted barcodes
    :type whitelist: str
    :param k: k-mer length
    :type k: int
    :return: Dictionary mapping all k-mers to
        indices in the whitelist corresponding to barcodes containing that k-mer
    :rtype: list, dict
    """
    kmer_to_bc_index = {}
    for index, bc in enumerate(whitelist):
        bc_kmers = find_kmers(bc, k)
        for bc_kmer in bc_kmers:
            if bc_kmer not in kmer_to_bc_index.keys():
                kmer_to_bc_index[bc_kmer] = set([index])
            else:
                kmer_to_bc_index[bc_kmer].add(index)

    return kmer_to_bc_index


def load_whitelist(whitelist_path, k):
    """
    Read in barcode whitelist and create dictionary mapping each k-mer to all
    barcodes in the whitelist containing that k-mer.

    :param whitelist_path: Path to the barcode whitelist
    :type whitelist_path: str
    :param k: k-mer length
    :type k: int
    :return: List of *sorted* whitelisted barcodes and dictionary mapping all k-mers to
        indices in the whitelist corresponding to barcodes containing that k-mer
    :rtype: list, dict
    """
    wl = []
    with open(whitelist_path) as file:
        for line in file:
            bc = line.strip().split("-")[0]
            wl.append(bc)

    wl.sort()

    return wl, generate_kmer_index(whitelist=wl, k=k)


def calc_leven_to_whitelist(bc_uncorr, whitelist, bc_len, kmer_to_bc_index):
    """
    Find minimum and runner-up barcode edit distance by iterating through the
    whitelist of expected barcodes.

    :param bc_uncorr: Uncorrected cell barcode
    :type bc_uncorr: str
    :param whitelist: Filtered whitelist of cell barcodes
    :type whitelist: list
    :param bc_len: length of barcode
    :type whitelist: int
    :param k: kmer length for whitelist indexing
    :type whitelist: int
    :return: Corrected barcode assignment, edit distance, and difference in edit
        distance between the top match and the next closest match
    :rtype: str, int, int
    """
    bc_corr = "-"
    bc_corr_leven = bc_len
    next_bc_corr_leven = bc_len

    # pull k value from kmer_to_bc_index
    k = len(list(kmer_to_bc_index.keys())[0])

    wl_filtered = filter_whitelist_by_kmers(
        wl=whitelist,
        kmers=find_kmers(bc_uncorr, k),
        kmer_to_bc_index=kmer_to_bc_index,
    )

    for wl_bc in wl_filtered:
        # d = ed.eval(bc_uncorr, wl_bc)  # Use the ed module here
        d = distance(bc_uncorr, wl_bc)  # levenshtein-python is much faster
        # d = hamming(bc_uncorr, wl_bc)

        if d < bc_corr_leven:
            next_bc_corr_leven = bc_corr_leven
            bc_corr_leven = d
            bc_corr = wl_bc
        elif d < next_bc_corr_leven:
            next_bc_corr_leven = d

        if d == 0:
            break

    next_match_diff = next_bc_corr_leven - bc_corr_leven

    return bc_corr, bc_corr_leven, next_match_diff


def process_tsv(
    tsv_in,
    tsv_out_full,
    tsv_out_slim,
    id_column,
    bc_columns,
    concat_bcs,
    whitelists,
    kmer_to_bc_indexes,
    max_levens,
    min_next_match_diffs,
    verbose=False,
    bc_update_counter=1000000,
):
    """
    Process a TSV file to correct barcodes based on a whitelist.

    Parameters:
    - tsv_in: Path to the input TSV file.
    - tsv_out_full: Path to the output TSV file with full information.
    - tsv_out_slim: Path to the output TSV file with slim information.
    - id_column: Column index for read IDs.
    - bc_columns: List of column indices for uncorrected barcodes.
    - concat_bcs: Boolean indicating whether to concatenate barcodes.
    - whitelists: Dictionary of whitelists for each barcode column.
    - kmer_to_bc_indexes: Dictionary of k-mer to barcode index mappings.
    - max_levens: List of maximum Levenshtein distances for correction.
    - min_next_match_diffs: List of minimum differences between first and second closest matches.
    - verbose: Boolean indicating whether to print verbose output.
    - bc_update_counter: Number of barcodes to process before printing an update.
    """
    null_bc_string = "-"

    # bc_update_counter=5000
    n_corrected = 0

    if verbose:
        processStartTime = time.time()
        print(f"Correcting barcode(s)..")
    # Prepare the output file
    with open(tsv_out_full, "w", newline="") as outfile_full:
        writer_full = csv.writer(outfile_full, delimiter="\t")

        with open(tsv_out_slim, "w", newline="") as outfile_slim:
            writer_slim = csv.writer(outfile_slim, delimiter="\t")

            read_count = 0
            with open(tsv_in, "r") as infile:
                reader = csv.reader(infile, delimiter="\t")

                for row in reader:
                    read_count += 1
                    if verbose and read_count % bc_update_counter == 0:
                        print(
                            f"{time.time()-processStartTime:.2f} - {read_count} reads processed..."
                        )

                    read_id = row[id_column]
                    barcodes = [row[i] for i in bc_columns]

                    row2write = [read_id]  # initialize row to write in output .tsv

                    if concat_bcs:
                        barcodes = ["".join(barcodes)]
                        RANGE = [
                            [
                                barcodes[0],
                                whitelists[0],
                                max_levens[0],
                                min_next_match_diffs[0],
                            ]
                        ]

                    else:
                        RANGE = [
                            [
                                barcodes[i],
                                whitelists[i],
                                max_levens[i],
                                min_next_match_diffs[i],
                            ]
                            for i in range(len(barcodes))
                        ]

                    # used for slim list for simple .bam tagging
                    concat_corrected_bc = []

                    for barcode, whitelist, max_leven, min_next_match_diff in RANGE:
                        st = time.time()

                        # skip reads w/ missing bc
                        # if barcode == null_bc_string:
                        #     continue

                        if len(barcode) != len(whitelist[0]):
                            print(
                                f"Error: Barcode length [{len(barcode)}] does not match whitelist [{len(whitelist[0])}] for read ID [{read_id}]"
                            )
                            sys.exit(1)

                        # quick check for perfect match
                        if barcode in whitelist:
                            corrected_bc = barcode
                            concat_corrected_bc.append(barcode)
                            bc_leven = 0
                            next_match_diff = "N"
                        else:
                            # calculate Levenshtein distance
                            (
                                corrected_bc,
                                bc_leven,
                                next_match_diff,
                            ) = calc_leven_to_whitelist(
                                bc_uncorr=barcode,
                                whitelist=whitelist,
                                bc_len=len(barcode),
                                kmer_to_bc_index=kmer_to_bc_indexes[i],
                            )
                            concat_corrected_bc.append(corrected_bc)

                            # k = len(list(kmer_to_bc_index.keys())[0])

                        # Write output(s)
                        if bc_leven == 0:
                            row2write.extend(
                                [barcode, corrected_bc, bc_leven, next_match_diff]
                            )
                        elif (
                            bc_leven <= max_leven
                            and next_match_diff >= min_next_match_diff
                        ):
                            row2write.extend(
                                [barcode, corrected_bc, bc_leven, next_match_diff]
                            )
                        else:
                            row2write.extend(
                                [barcode, null_bc_string, bc_leven, next_match_diff]
                            )

                    # end row-wise BC correction
                    if null_bc_string in concat_corrected_bc:
                        # only write fully correctable bcs to slim output
                        writer_full.writerow(row2write)
                    else:
                        n_corrected += 1
                        writer_full.writerow(row2write)
                        writer_slim.writerow([read_id, "".join(concat_corrected_bc)])
    if n_corrected == 0:
        print(f"WARNING - NO BARCODES MATCHED THE WHITELIST")
    if verbose:
        print(f"{time.time()-processStartTime:.2f} - Finished correcting!")


if __name__ == "__main__":
    args = parse_args()

    # Check if input TSV file exists
    if not os.path.isfile(args.tsv_in):
        print(f"Error: Input TSV file '{args.tsv_in}' does not exist.")
        sys.exit(1)

    # Check if whitelist files exist
    for whitelist_file in args.whitelist_files:
        if not os.path.isfile(whitelist_file):
            print(f"Error: Whitelist file '{whitelist_file}' does not exist.")
            sys.exit(1)

    # param checks ------
    if len(args.max_levens) != len(args.bc_columns) and len(args.max_levens) == 1:
        args.max_levens = rep(val=args.max_levens[0], n=len(args.bc_columns))

    if (
        len(args.min_next_match_diffs) != len(args.bc_columns)
        and len(args.min_next_match_diffs) == 1
    ):
        args.min_next_match_diffs = rep(
            val=args.min_next_match_diffs[0], n=len(args.bc_columns)
        )

    # Add check for whitelist files matching bc_columns when concat_bcs is False
    if not args.concat_bcs and len(args.whitelist_files) != len(args.bc_columns):
        print(f"Error: Number of whitelist files ({len(args.whitelist_files)}) does not match the number of barcode columns ({len(args.bc_columns)}) when concat_bcs is set to False.")
        sys.exit(1)

    # Check if column indices are valid
    with open(args.tsv_in, "r") as infile:
        reader = csv.reader(infile, delimiter="\t")
        header = next(reader)
        num_columns = len(header)

        if args.id_column >= num_columns:
            print(f"Error: Read ID column index ({args.id_column}) is out of range.")
            sys.exit(1)

        for bc_column in args.bc_columns:
            if bc_column >= num_columns:
                print(f"Error: Barcode column index ({bc_column}) is out of range.")
                sys.exit(1)

    # Print run settings for log files ----
    print(
        f"Input tsv:                        {args.tsv_in}\n"
        f"Output tsv (Full Info):           {args.tsv_out_full}\n"
        f"Output tsv (CBs only):            {args.tsv_out_slim}\n"
        f"Read ID column:                   {args.id_column}\n"
        f"Barcode column(s):                {args.bc_columns}\n"
        f"Concatenate uncorrected BCs?:     {args.concat_bcs}\n"
        f"Whitelist file(s):                {args.whitelist_files}\n"
        f"Whitelist kmer filter length(k):  {args.k}\n"
        f"Maximum levenshtein distance(s):  {args.max_levens}\n"
        f"Second match difference(s):       {args.min_next_match_diffs}\n"
        f"Number of threads:                {args.threads}\n"
    )

    if args.concat_bcs and len(list(args.whitelist_files)) > 1:
        print(f"Need a merged barcode whitelist!")
        sys.exit(1)

    # prep whitelists
    if len(args.whitelist_files) != len(args.bc_columns):
        print(f"Using this whitelist for all corrections:   {args.whitelist_files[0]}")
        args.whitelist_files = rep(args.whitelist_files[0], len(args.bc_columns))

    whitelists = {}
    kmer_to_bc_indexes = {}
    for i, whitelist_file in enumerate(args.whitelist_files):
        wl, kmer_to_bc_index = load_whitelist(whitelist_file, k=args.k)
        # TODO- auto set k based on bc length (dependent on seq error rate)

        whitelists[i] = wl
        kmer_to_bc_indexes[i] = kmer_to_bc_index

    # Single-threaded = verbose
    if args.threads == 1:
        print(f"{currentTime()} - Running on {args.threads} thread...")
        print("")

        process_tsv(
            tsv_in=args.tsv_in,
            tsv_out_full=args.tsv_out_full,
            tsv_out_slim=args.tsv_out_slim,
            id_column=args.id_column,
            bc_columns=args.bc_columns,
            concat_bcs=args.concat_bcs,
            whitelists=whitelists,
            kmer_to_bc_indexes=kmer_to_bc_indexes,
            max_levens=args.max_levens,
            min_next_match_diffs=args.min_next_match_diffs,
            verbose=True,
        )

    # Multi-threaded = not verbose
    elif args.threads > 1:
        print(f"{currentTime()} - Running on {args.threads} threads...")
        # Source: https://superfastpython.com/multiprocessing-pool-for-loop/
        import multiprocessing

        print("")
        print(f"{currentTime()} - Splitting input tsv...")
        temp_dir = generate_temp_dir_name()

        if not os.path.exists(temp_dir):
            os.makedirs(temp_dir)

        print(f"Temporary directory: {temp_dir}")
        temp_tsvs_in, n_bcs = split_file(
            file_path=args.tsv_in, temp_dir=temp_dir, num_chunks=args.threads
        )

        print(f"{currentTime()} - Correcting {n_bcs} barcodes...")

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
                args.id_column,
                args.bc_columns,
                args.concat_bcs,
                whitelists,
                kmer_to_bc_indexes,
                args.max_levens,
                args.min_next_match_diffs,
            )
            for i in list(range(args.threads))
        ]

        with multiprocessing.Pool(args.threads) as pool:
            multi_out = pool.starmap(process_tsv, items)

        print("")
        print(f"{currentTime()} - Concatenating results...")
        concatenate_files(
            filenames=temp_tsvs_out_full, output_filename=args.tsv_out_full
        )
        concatenate_files(
            filenames=temp_tsvs_out_slim, output_filename=args.tsv_out_slim
        )

        print(f"{currentTime()} - Removing temp files (temp_dir: `{temp_dir}`)...")
        os.system(f"rm -rf {temp_dir}") if os.path.exists(temp_dir) else None

        print(f"{currentTime()} - Done!")
        print("")
        print(f"Finished in {time.time() - startTime:.2f} seconds")
    else:
        print(f"Incorrect number of threads specified [{args.threads}]...")
