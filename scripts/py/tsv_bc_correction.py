import sys
import csv
import argparse
from itertools import product
import editdistance as ed
import sys
import csv
import argparse
from itertools import product


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


def main(args):
    input_file = args.input_file
    whitelist_files = args.whitelist_files.split(",")
    min_dist = int(args.min_dist) if args.min_dist else 2

    # Prepare the output file
    with open(args.output_file, "w", newline="") as outfile:
        writer = csv.writer(outfile, delimiter="\t")
        # writer.writerow(['Read_ID', 'Original_Barcode', 'Corrected_Barcode', 'Hamming_Distance'])

        with open(input_file, "r") as infile:
            reader = csv.reader(infile, delimiter="\t")
            next(reader)  # Skip header row

            for row in reader:
                read_id = row[0]
                barcodes = row[1:]

                whitelists = []
                for i, barcode in enumerate(barcodes):
                    whitelist_file = whitelist_files[i]
                    wl, kmer_to_bc_index = load_whitelist(
                        whitelist_file, k=5
                    )  # Assuming k=5 for simplicity
                    whitelists.append(wl)

                row2write = [read_id]
                for barcode, whitelist in zip(barcodes, whitelists):
                    corrected_bc, bc_match_ed, next_match_diff = calc_ed_with_whitelist(
                        barcode, whitelist
                    )
                    if bc_match_ed < min_dist:
                        # row2write.append(f"{barcode}\t{corrected_bc}\t{bc_match_ed}")
                        row2write.extend(
                            [barcode, corrected_bc, bc_match_ed, next_match_diff]
                        )
                    else:
                        # row2write.append(f"{barcode}\t \t{bc_match_ed}")
                        row2write.extend([barcode, " ", bc_match_ed, next_match_diff])
                writer.writerow(row2write)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Correct barcodes in a TSV file.")
    parser.add_argument(
        "--input_file", required=True, help="Path to the input TSV file."
    )
    parser.add_argument(
        "--output_file",
        required=True,
        help="Path to the output TSV file. Columns contain ['Read_ID', 'Original_Barcode', 'Corrected_Barcode', 'Hamming_Distance']",
    )
    parser.add_argument(
        "--whitelist_files",
        required=True,
        help="Comma-separated list of whitelist file paths.",
    )
    parser.add_argument(
        "--min_dist",
        type=int,
        default=2,
        help="Minimum Hamming distance for correction (default: 2).",
    )

    args = parser.parse_args()
    main(args)
