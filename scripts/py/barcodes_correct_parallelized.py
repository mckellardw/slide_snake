import argparse
import collections
import gzip
import logging
import math
import multiprocessing
import os
import pathlib
import re
import shutil
import tempfile

import editdistance as ed

# import parasail
import pysam
from tqdm import tqdm

logger = logging.getLogger(__name__)

# for debugging
verbose = False  # True


def parse_args():
    # Create argument parser
    parser = argparse.ArgumentParser()

    # Positional mandatory arguments
    parser.add_argument(
        "bam",
        help="Sorted BAM file of stranded sequencing reads aligned to a reference. \
            Alignments must have the CR and CY tags.",
        type=str,
    )

    parser.add_argument(
        "whitelist", help="File containing list of expected cell barcodes", type=str
    )

    # Optional arguments
    parser.add_argument(
        "-t", "--threads", help="Threads to use [4]", type=int, default=4
    )

    parser.add_argument(
        "-k", help="Kmer size to use for whitelist filtering [5]", type=int, default=5
    )

    parser.add_argument(
        "--output_bam",
        help="Output BAM file containing aligned reads with tags for uncorrected \
        barcodes (CR), corrected barcodes (CB), barcode QVs (CY), uncorrected \
        UMIs (UR), and UMI QVs (UY) [bc_corr.umi_uncorr.sorted.bam]",
        type=str,
        default="bc_corr.umi_uncorr.sorted.bam",
    )

    parser.add_argument(
        "--output_counts",
        help="Output TSV file containing counts for each of the assigned \
        barcodes [barcode_counts.tsv]",
        type=str,
        default="barcode_counts.tsv",
    )

    parser.add_argument(
        "--max_ed",
        help="Max edit distance between putative barcode \
                        and the matching whitelist barcode [2]",
        type=int,
        default=2,
    )

    parser.add_argument(
        "--min_ed_diff",
        help="Min difference in edit distance between the \
                        (1) putative barcode vs top hit and (2) putative \
                        barcode vs runner-up hit [2]",
        type=int,
        default=2,
    )

    parser.add_argument(
        "--kit",
        help="Specify either the 10X 3' gene expression kit (3prime), the 5' \
        gene expression kit (5prime), or the multiome kit (multiome) This \
        determines which adapter sequences to search for in the reads \
        [3prime]",
        default="3prime",
        type=str,
    )

    parser.add_argument(
        "--adapter1_suff_length",
        help="Use this many suffix bases from adapter1 sequence  in the \
            alignment query. For example, specifying 12 would mean that the last \
            12 bases of the specified read1 sequence will be included in the \
            probe sequence [10]",
        default=10,
        type=int,
    )

    parser.add_argument(
        "-T",
        "--polyT_length",
        help="Length of polyT sequence to use in the alignment query (ignored \
        with --kit=5prime) [10]",
        type=int,
        default=10,
    )

    # parser.add_argument(
    #     "--barcode_length", help="Cell barcode length [16]", type=int, default=16
    # )

    parser.add_argument("--umi_length", help="UMI length [12]", type=int, default=12)

    parser.add_argument(
        "-w",
        "--window",
        help="Number of bases to query at start of read [100]",
        type=int,
        default=100,
    )

    parser.add_argument(
        "-o", "--gap_open", help="Gap open penalty [2]", type=int, default=2
    )

    parser.add_argument(
        "-e", "--gap_extend", help="Gap extend penalty [4]", type=int, default=4
    )

    parser.add_argument("-m", "--match", help="Match score [5]", type=int, default=5)

    parser.add_argument(
        "-x", "--mismatch", help="Mismatch score [-1]", type=int, default=-1
    )

    parser.add_argument(
        "-n",
        "--acg_to_n_match",
        help="Score for A/C/G<-->N match [1]",
        type=int,
        default=1,
    )

    parser.add_argument(
        "-s", "--t_to_n_match", help="Score for T<-->N match [1]", type=int, default=1
    )

    parser.add_argument(
        "--verbosity",
        help="logging level: <=2 logs info, <=3 logs warnings",
        type=int,
        default=2,
    )

    # Parse arguments
    args = parser.parse_args()

    # verify kit selection
    if (args.kit != "3prime") and (args.kit != "5prime"):
        raise Exception("Invalid kit name! Specify either '3prime' or '5prime'.")

    if args.kit == "3prime":
        # Read1 adapter
        args.adapter1_seq = "CTACACGACGCTCTTCCGATCT"
        # TSO adapter
        args.adapter2_seq = "ATGTACTCTGCGTTGATACCACTGCTT"
    elif args.kit == "5prime":
        # Read1 adapter
        args.adapter1_seq = "CTACACGACGCTCTTCCGATCT"
        # Poly-dT RT adapter
        args.adapter2_seq = "GTACTCTGCGTTGATACCACTGCTT"

    # Create temp dir and add that to the args object
    p = pathlib.Path(args.output_bam)
    tempdir = tempfile.TemporaryDirectory(prefix="tmp.", dir=p.parents[0])
    args.tempdir = tempdir.name

    return args


def init_logger(args):
    logging.basicConfig(
        format="%(asctime)s -- %(message)s", datefmt="%Y-%m-%d %H:%M:%S"
    )
    logging_level = args.verbosity * 10
    logging.root.setLevel(logging_level)
    logging.root.handlers[0].addFilter(lambda x: "NumExpr" not in x.msg)


# TODO potentially could replace with [fast-edit-distance](https://github.com/youyupei/fast_edit_distance)
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
        d = ed.eval(bc_uncorr, wl_bc)
        if d < bc_match_ed:
            next_bc_match_ed = bc_match_ed
            bc_match_ed = d
            bc_match = wl_bc
        elif d < next_bc_match_ed:
            next_bc_match_ed = d
    next_match_diff = next_bc_match_ed - bc_match_ed

    return bc_match, bc_match_ed, next_match_diff


# Removed `chrom` input
def process_bam_records(tup):
    """
    Process BAM records to assign each read a corrected cell barcode and an
    uncorrected UMI. Do this by loading and processing the barcode whitelist
    then iterating over alignments.
    For each alignment:
    1. Calculate edit distance between uncorrected barcode and barcodes in whitelist
    2.

    :param tup: Tuple containing the input arguments
    :type tup: tup
    :return: Path to a temporary BAM file
    :rtype: str
    """
    input_bam = tup[0]
    output_bam = tup[1]
    threads = tup[2]
    args = tup[3]

    verbose = False  # for debugging

    # Load barcode whitelist and map kmers to indices in whitelist for faster barcode matching
    whitelist, kmer_to_bc_index = load_whitelist(args.whitelist, args.k)

    # Open input BAM file
    bam = pysam.AlignmentFile(input_bam, "rb")

    # Write temp file or straight to output file depending on use case
    if threads > 1:
        # Open temporary output BAM file for writing
        suff = os.path.basename(input_bam)
        chrom_bam = tempfile.NamedTemporaryFile(
            prefix="tmp.align.", suffix=suff, dir=args.tempdir, delete=False
        )
        bam_out_fn = chrom_bam.name
    else:
        bam_out_fn = output_bam

    bam_out = pysam.AlignmentFile(output_bam, "wb", template=bam)

    barcode_counter = collections.Counter()

    for align in bam.fetch():  # contig=chrom
        # Make sure each alignment in this BAM has an uncorrected barcode and barcode QV
        # TODO add quality score extraction in barcode/UMI extraction step
        assert align.has_tag("CR"), "CR tag not found"  # and align.has_tag("CY")

        bc_uncorr = align.get_tag("CR")

        # Don't consider any uncorrected barcodes that are shorter than k
        if len(bc_uncorr) >= args.k:
            # Decompose uncorrected barcode into N k-mers
            bc_uncorr_kmers = split_seq_into_kmers(bc_uncorr, args.k)

            # Filter the whitelist to only those with at least one of the k-mers
            # from the uncorrected barcode
            filt_whitelist = filter_whitelist_by_kmers(
                whitelist, bc_uncorr_kmers, kmer_to_bc_index
            )

            # Calc edit distances between uncorrected barcode and the filtered
            # whitelist barcodes
            bc_match, bc_match_ed, next_match_diff = calc_ed_with_whitelist(
                bc_uncorr, filt_whitelist
            )

            # For debugging/optimizing:
            # if verbose:
            #     print(f"bc_match = {bc_match}")
            #     print(f"bc_match_ed = {bc_match_ed}")
            #     print(f"next_match_diff = {next_match_diff}")

            # Check barcode match edit distance and difference to runner-up edit distance
            condition1 = bc_match_ed <= args.max_ed
            condition2 = next_match_diff >= args.min_ed_diff
            if condition1 and condition2:
                # Add corrected cell barcode = CB:Z
                align.set_tag("CB", bc_match, value_type="Z")

            # Only write BAM entry in output file if we've assigned a corrected
            # barcode and an uncorrected UMI
            if align.has_tag("CB") and align.has_tag("UR"):
                bam_out.write(align)
                barcode_counter[bc_match] += 1

    bam.close()
    bam_out.close()

    return bam_out_fn, barcode_counter


def launch_pool(func, func_args, procs=1):
    """
    Use multiprocessing library to create pool and map function calls to
    that pool

    :param procs: Number of processes to use for pool
    :type procs: int, optional
    :param func: Function to exececute in the pool
    :type func: function
    :param func_args: List containing arguments for each call to function <funct>
    :type func_args: list
    :return: List of results returned by each call to function <funct>
    :rtype: list
    """
    p = multiprocessing.Pool(processes=procs)
    try:
        results = list(tqdm(p.imap(func, func_args), total=len(func_args)))
        p.close()
        p.join()
    except KeyboardInterrupt:
        p.terminate()
    return results


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


def get_bam_info(bam):
    """
    Use `samtools idxstat` to get number of alignments and names of all contigs
    in the reference.

    :param bam: Path to sorted BAM file
    :type bame: str
    :return: Sum of all alignments in the BAM index file and list of all chroms
    :rtype: int,list
    """
    bam = pysam.AlignmentFile(bam, "rb")
    stats = bam.get_index_statistics()
    n_aligns = int(sum([contig.mapped for contig in stats]))  # number of aligned reads
    n_unmapped = bam.nocoordinate  # number of UNaligned reads
    n_reads = n_aligns + n_unmapped  # number of total reads
    chroms = dict(
        [(contig.contig, contig.mapped) for contig in stats if contig.mapped > 0]
    )
    bam.close()
    return n_reads, n_aligns, chroms


# DWM
def split_bam_file(input_bam, num_files, n_reads=None, output_dir=None, index=True):
    """
    Splits a BAM file into N files and returns the names of the split BAM files.

    Parameters:
    - input_bam: Path to the input BAM file.
    - num_files: Number of files to split the input BAM into.
    - n_reads: Optional. Number of reads in the BAM file. If not provided, the function will count the reads.
    - output_prefix: Optional. Prefix for the output BAM files. If not provided, the function will use the directory of the input BAM file.
    - indexL Boolean. Whether or not to index the split bam files. Default: True

    Returns:
    - A list of the names of the split BAM files.
    """
    # Handle the output directory, make it if it doesn't exist
    if output_dir is None:
        output_dir = os.path.splitext(input_bam)[0]
    else:
        output_dir = os.path.dirname(output_dir)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Open the input BAM file
    with pysam.AlignmentFile(input_bam, "rb") as input_bam_file:
        # If n_reads is not provided, count the reads
        if n_reads is None:
            n_reads = input_bam_file.count()  # this does not include unaligned

        # Calculate the number of reads per file
        reads_per_file = math.ceil(n_reads / num_files)

        # Initialize counters and file handles
        current_read_count = 0
        current_file_number = 1
        current_output_bam = pysam.AlignmentFile(
            f"{output_dir}/{current_file_number:03d}.bam", "wb", template=input_bam_file
        )
        split_files = []  # List to store the names of the split BAM files

        # Iterate over the reads in the input BAM file
        for read in input_bam_file:
            # If we've reached the desired number of reads for the current file,
            #    close the current output file and open a new one
            if current_read_count >= reads_per_file:
                if current_output_bam is not None:
                    current_output_bam.close()
                current_output_bam = pysam.AlignmentFile(
                    f"{output_dir}/{current_file_number:03d}.bam",
                    "wb",
                    template=input_bam_file,
                )
                split_files.append(
                    f"{output_dir}/{current_file_number:03d}.bam"
                )  # Add the current output file name to the list
                current_file_number += 1
                current_read_count = 0

            # Write the current read to the current output file
            current_output_bam.write(read)
            current_read_count += 1

        # Close the last output file
        if current_output_bam is not None:
            current_output_bam.close()
            # split_files.append(f"{output_dir}/{current_file_number:03d}.bam") # Add the last output file name to the list

    if index:
        logger.info("Indexing split BAMs...")
        for bam in split_files:
            pysam.index(bam)

    return split_files


def main(args):
    init_logger(args)
    logger.info(f"Using {args.threads} thread(s)")
    logger.info("Getting BAM statistics")
    n_reads, n_aligns, chroms = get_bam_info(args.bam)
    logger.info(f"    Found {n_reads} reads in {args.bam}")
    logger.info(f"    Found {n_aligns} alignments in {args.bam}")

    logger.info(f"Correcting barcodes in {args.bam}")

    if args.threads > 1:
        # Create temporary directory
        if os.path.exists(args.tempdir):
            shutil.rmtree(args.tempdir, ignore_errors=True)
        os.mkdir(args.tempdir)

        # Split bam and process across multiple threads
        logger.info(f"    Splitting BAM into {args.threads} BAMs...")

        split_bam_files = split_bam_file(
            input_bam=args.bam,
            # output_prefix=,
            num_files=args.threads,  # 1 bam per thread
            n_reads=n_reads,
        )

        func_args = []
        for bam in split_bam_files:
            func_args.append((bam, bam.replace(".bam", "_corrected.bam"), 1, args))
            logger.info(f"    {bam}")

        # if verbose: # for debuging
        #     logger.info(f"func_args for pool run:")
        #     for farg in func_args:
        #         logger.info(f"    {farg}")

        results = launch_pool(
            func=process_bam_records, func_args=func_args, procs=args.threads
        )
        chrom_bam_fns, barcode_counters = zip(*results)

        barcode_counter = sum(barcode_counters, collections.Counter())

        tmp_bam = tempfile.NamedTemporaryFile(
            prefix="tmp.align.", suffix=".unsorted.bam", dir=args.tempdir, delete=False
        )
        merge_parameters = ["-f", tmp_bam.name] + list(chrom_bam_fns)
        pysam.merge(*merge_parameters)

        pysam.sort("-@", str(args.threads), "-o", args.output_bam, tmp_bam.name)

        logger.info("Cleaning up temporary files")
        shutil.rmtree(args.tempdir, ignore_errors=True)
        shutil.rmtree(  # remove split bams
            os.path.dirname(split_bam_files[0]), ignore_errors=True
        )
    else:
        func_args = (args.bam, args)  # chrom,
        chrom_bam_fn, barcode_counter = process_bam_records(func_args)
        logger.info(f"Written to {chrom_bam_fn}")

    with open(args.output_counts, "w") as f:
        for bc, n in barcode_counter.most_common():
            f.write(f"{bc}\t{n}\n")

    logger.info("Done.")


if __name__ == "__main__":
    args = parse_args()

    # Find the barcode length from the whitelist instead of requiring the argument
    ## Not currently used in script?
    args.barcode_length = len(open(args.whitelist).readline())

    main(args)
