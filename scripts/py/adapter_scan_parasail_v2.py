import argparse
import gzip
import multiprocessing
import os
import pathlib
import shutil
import sys
import tempfile
import time
import io

import numpy as np
import pandas as pd
import parasail
import pysam


# Example usage:
"""
python scripts/py/adapter_scan_parasail.py \
    --fq_in "sandbox/adapter_scan/maxima/merged.fq.gz" \
    --fq_out "sandbox/adapter_scan/maxima/merged_stranded.fq.gz" \
    --output_tsv "sandbox/adapter_scan/maxima/adapter_scan.tsv" \
    --threads 56 \
    --batch_size 100000 \
    --adapter1_seq "CTACACGACGCTCTTCCGATCT" \
    --adapter2_seq "ATGTACTCTGCGTTGATACCACTGCTT" \
    --min_adapter_match 0.7 
"""


def currentTime():
    """Return the current time formatted as 'MM/DD/YY | HH:MM:SS'."""
    return time.strftime("%D | %H:%M:%S", time.localtime())


def print_log(msg):
    """Print a log message with the current time."""
    print(f"{currentTime()} | {msg}")


def print_error(msg):
    """Print an error message with the current time."""
    print(f"{currentTime()} | ERROR: {msg}", file=sys.stderr)
    sys.exit(1)


def parse_args():
    """Parse command-line arguments and create a temporary directory."""
    parser = argparse.ArgumentParser()

    # Mandatory arguments
    parser.add_argument("--fq_in", help="FASTQ of ONT reads", type=str)

    # Optional arguments
    parser.add_argument(
        "--fq_out",
        help="Output file name for (gzipped) stranded FASTQ entries [stranded.fastq.gz]",
        type=str,
        default="stranded.fastq.gz",
    )

    parser.add_argument(
        "--output_tsv",
        help="Output file name for adapter configurations [adapters.tsv]",
        type=str,
        default="adapters.tsv",
    )

    parser.add_argument(
        "-t", "--threads", help="Threads to use [4]", type=int, default=4
    )

    parser.add_argument(
        "-b",
        "--batch_size",
        help="Number of reads per batch [100000]",
        type=int,
        default=100000,
    )

    parser.add_argument(
        "-f",
        "--adapter1_seq",
        help="Forward primer sequences (Read1)- this is usually followed by a poly(T) region, \
            reverse-complement of the capture sequence. Multiple sequences can be passed for more \
            rigorous stranding, or if iso-primers are used. If you use multiple sequences, the \
            order should reflect their 5'-3' order on the read. \
            For TXG, Curio, etc. use [CTACACGACGCTCTTCCGATCT]",
        nargs="+",
        default=["CTACACGACGCTCTTCCGATCT"],
        type=str,
    )

    parser.add_argument(
        "-r",
        "--adapter2_seq",
        help="Reverse complement of the primer sequences (rcTSO, minus the CCC prefix). \
            For TXG/Curio chemistries, use [ATGTACTCTGCGTTGATACCACTGCTT] (default). \
            For uMRT chemistries, use [GAGAGAGGAAAGAGAGAGAGAGGG]",
        nargs="+",
        default=None,
        type=str,
    )

    parser.add_argument(
        "-i",
        "--min_adapter_match",
        help="Minimum adapter alignment identity for parasail [0.7]",
        type=float,
        default=0.7,
    )

    # Parse arguments
    args = parser.parse_args()

    # If adapter2_seq was not provided, set to empty list and enable single adapter mode
    if not args.adapter2_seq:
        args.adapter2_seq = []
        args.single_end_mode = True
    else:
        args.single_end_mode = False

    # Create temp dir and add that to the args object
    p = pathlib.Path(args.output_tsv)
    tempdir = tempfile.TemporaryDirectory(prefix="tmp.", dir=p.parents[0])
    args.tempdir = tempdir.name

    return args


# complement translation table with support for regex punctuation
COMPLEMENT_TRANS = str.maketrans(
    "ACGTWSMKRYBDHVNacgtwsmkrybdhvn", "TGCAWSKMYRVHDBNtgcawskmyrvhdbn"
)


# TODO- allow support for multi-alignments (important for TSO/PCR concatemers!)
def align_parasail(read, adapter, min_adapter_match=0.7, matrix=None):
    """
    Align adapter sequence to a read using parasail (Smith-Waterman local alignment)
    - source: https://github.com/jeffdaily/parasail-python

    Parameters:
        read (str): The DNA sequence of the read to which the adapter will be aligned.
        adapter (str): The DNA sequence of the adapter to align to the read.
        min_adapter_match (float): The minimum adapter identity threshold (between 0 and 1)
                                to consider the alignment valid.

    Returns:
        parasail.Alignment or None: The alignment object if the alignment score meets
                                     the threshold, otherwise None.
    """
    # Create scoring matrix and perform alignment
    if matrix is None:
        matrix = parasail.matrix_create("ACGT", 1, -1)

    # s1=query, s2=target
    alignment = parasail.sw_trace_striped_32(
        s1=adapter, s2=read, open=10, extend=1, matrix=matrix
    )

    # Calculate mismatches and check if alignment score meets threshold
    min_score = round(min_adapter_match * len(adapter))
    if alignment.score >= min_score:
        return alignment
    return None


def add_entry_to_read_info(
    read_info,
    orig_read_id,
    new_read_id,
    readlen,
    start,
    end,
    strand,
    adapter_config,
    lab,
):
    """Add processed read information to the read_info dictionary."""
    read_info[orig_read_id][new_read_id] = {
        "readlen": readlen,
        "start": start,
        "end": end,
        "strand": strand,
        "lab": lab,
        "adapter_config": adapter_config,
    }
    return read_info


def determine_strand_from_parasail(
    align_df,
    single_end_mode,
    read_seq,
    read_id,
    read_info,
):
    """
    Determine the strand orientation and update the read_info dictionary.

    Args:
        align_df (pd.DataFrame): Dataframe containing alignment information.
        single_end_mode (bool): Whether single adapter mode is enabled.
        read_seq (str): The read sequence.
        read_id (str): The read identifier.
        read_info (dict): Dictionary to store processed read information.

    Returns:
        dict: Updated read_info dictionary.

    Configuration options:
        "full_len",
        "adapter1_single",
        "adapter1_double",
        "adapter2_single",
        "adapter2_double",
        "other",
        "other_ambiguous",
        "no_adapters"
    """
    # Sort alignments by position
    align_df = align_df.sort_values("qilo").reset_index(drop=True)

    # get strand(s)
    strands = align_df["qstrand"].tolist()

    # get adapter_config
    adapter_config = "-".join(align_df["query"].tolist())
    new_read_id = align_df["target"].iloc[0]

    if single_end_mode:
        # single adapter mode
        if all(s == "+" for s in strands):
            strand = "+"
            start = min(align_df["qihi"].tolist())
            end = max(align_df["qilo"].tolist())
            readlen = end - start
            lab = "full_len"
        elif all(s == "-" for s in strands):
            strand = "-"
            start = min(align_df["qihi"].tolist())
            end = len(read_seq)
            readlen = end - start
            lab = "adapter1_single"
        else:
            strand = "="
    else:  # normal paired adapters mode
        if all(s == "+" for s in strands):
            strand = "+"
            match adapter_config:
                case "adapter1_f_0-adapter2_r_0":
                    start = min(align_df["qihi"].tolist())
                    end = max(align_df["qilo"].tolist())
                    readlen = end - start
                    lab = "full_len"
                case "adapter1_f_0-adapter1_r_0":
                    start = min(align_df["qihi"].tolist())
                    end = max(align_df["qilo"].tolist())
                    readlen = end - start
                    lab = "double_adapter1"
                case "adapter1_f_0":
                    start = min(align_df["qihi"].tolist())
                    end = len(read_seq)
                    readlen = end - start
                    lab = "adapter1_single"
                case "adapter2_r_0":
                    start = 0
                    end = max(align_df["qilo"].tolist())
                    readlen = end - start
                    lab = "adapter2_single"
                case _:
                    start = 0
                    end = len(read_seq)
                    readlen = len(read_seq)
                    lab = "other"
        elif all(s == "-" for s in strands):  # TODO- flip start/end for RCed reads!
            strand = "-"
            match adapter_config:
                case "adapter2_f_0-adapter1_r_0":
                    # simple neg case
                    start = min(align_df["qihi"].tolist())
                    end = max(align_df["qilo"].tolist())
                    readlen = end - start
                    lab = "full_len"
                case "adapter1_r_0":
                    start = min(align_df["qihi"].tolist())
                    end = len(read_seq)
                    readlen = end - start
                    lab = "adapter1_single"
                case "adapter2_f_0":
                    start = 0
                    end = max(align_df["qilo"].tolist())
                    readlen = end - start
                    lab = "adapter2_single"
                case _:
                    start = 0
                    end = len(read_seq)
                    readlen = len(read_seq)
                    lab = "other"
        else:
            # get scores for each strand
            pos_scores = np.sum(align_df[align_df["qstrand"] == "+"]["score"].tolist())
            neg_scores = np.sum(align_df[align_df["qstrand"] == "-"]["score"].tolist())

            ratio_threshold = 1.1  # summed alignment score ratio for [confidently] determining strand- arbitrary...

            if pos_scores > ratio_threshold * neg_scores:
                strand = "+"
            elif neg_scores > ratio_threshold * pos_scores:
                strand = "-"
            else:
                strand = "="  # ambiguous strand

            match adapter_config:
                case "adapter1_f_0-adapter1_r_0" | "adapter1_r_0-adapter1_f_0":
                    start = min(align_df["qihi"].tolist())
                    end = max(align_df["qilo"].tolist())
                    readlen = end - start
                    lab = "adapter1_double"
                case "adapter2_f_0-adapter2_r_0" | "adapter2_r_0-adapter2_f_0":
                    start = min(align_df["qihi"].tolist())
                    end = max(align_df["qilo"].tolist())
                    readlen = end - start
                    lab = "adapter2_double"
                # TODO- cases for concatemers - eg: "adapter1_f_0-adapter2_r_0-adapter2_f_0-adapter1_r_0"
                case _:
                    start = 0
                    end = len(read_seq)
                    readlen = len(read_seq)
                    lab = "other_ambiguous"

        read_info = add_entry_to_read_info(
            read_info=read_info,
            orig_read_id=read_id,
            new_read_id=new_read_id,
            readlen=readlen,
            start=start,
            end=end,
            strand=strand,
            adapter_config=adapter_config,
            lab=lab,
        )

    return read_info


def parse_parasail(
    reads,
    single_end_mode,
    adapter1_seq,
    adapter2_seq,
    min_adapter_match,
):
    """Process reads to find adapter alignments and valid adapter pairs; return the read_info dict."""
    read_info = {}
    for read in reads:
        read_id = read.name
        read_seq = read.sequence
        read_info[read_id] = {}
        alignments = []

        # Align adapter1 sequences
        for i, adapter in enumerate(adapter1_seq):
            # fwd alignment
            alignment = align_parasail(
                read=read_seq, adapter=adapter, min_adapter_match=min_adapter_match
            )
            if alignment:
                alignments.append(
                    {
                        "target": read_id,
                        "query": f"adapter1_f_{i}",
                        "qilo": alignment.end_ref - len(adapter) + 1,  # start
                        "qihi": alignment.end_ref + 1,  # end
                        "qstrand": "+",
                        "score": alignment.score,
                        "query_seq": adapter,
                    }
                )

            # Reverse complement alignment
            alignment_rc = align_parasail(
                read_seq,
                adapter[::-1].translate(COMPLEMENT_TRANS),
                min_adapter_match,
            )
            if alignment_rc:
                alignments.append(
                    {
                        "target": read_id,
                        "query": f"adapter1_r_{i}",
                        "qilo": alignment_rc.end_ref - len(adapter) + 1,
                        "qihi": alignment_rc.end_ref + 1,
                        "qstrand": "-",
                        "score": alignment_rc.score,
                        "query_seq": adapter[::-1].translate(COMPLEMENT_TRANS),
                    }
                )

        # Align adapter2 sequences if not in single adapter mode
        if not single_end_mode:
            for i, adapter in enumerate(adapter2_seq):
                # fwd alignment
                result = align_parasail(
                    read=read_seq, adapter=adapter, min_adapter_match=min_adapter_match
                )
                if result:
                    alignments.append(
                        {
                            "target": read_id,
                            "query": f"adapter2_f_{i}",
                            "qilo": result.end_ref - len(adapter) + 1,
                            "qihi": result.end_ref + 1,
                            "qstrand": "-",
                            "score": result.score,
                            "query_seq": adapter,
                        }
                    )

                # Reverse complement alignment
                alignment_rc = align_parasail(
                    read_seq,
                    adapter[::-1].translate(COMPLEMENT_TRANS),
                    min_adapter_match,
                )
                if alignment_rc:
                    alignments.append(
                        {
                            "target": read_id,
                            "query": f"adapter2_r_{i}",
                            "qilo": alignment_rc.end_ref - len(adapter) + 1,
                            "qihi": alignment_rc.end_ref + 1,
                            "qstrand": "+",
                            "score": alignment_rc.score,
                            "query_seq": adapter[::-1].translate(COMPLEMENT_TRANS),
                        }
                    )
        # Process alignments to determine valid adapter pairs
        if alignments:
            align_df = pd.DataFrame(alignments)  # one alignment per row

            #### for debugging
            # save_dict = {
            #     "41061408-bbf4-4d9e-9835-a01a7fa3bee7": "full_len",
            #     "76364f8f-ac7a-4773-aeac-049ee4fbdc6e": "other",
            #     "14e386f9-c658-45f9-9ace-760fe5f9b2e2": "ambiguous",
            #     "659ce72a-7fb2-4b1c-b71a-b471ba9c0b7d": "other",
            # }
            # if read_id in save_dict.keys():
            #     output_dir = f"sandbox/adapter_scan/alignments/{save_dict[read_id]}"
            #     os.makedirs(output_dir, exist_ok=True)
            #     align_df.to_csv(
            #         f"{output_dir}/{read_id}.tsv",
            #         sep="\t",
            #         index=False,
            #     )

            read_info = determine_strand_from_parasail(
                align_df,
                single_end_mode,
                read_seq,
                read_id,
                read_info,
            )
        else:
            # add ambiguous entry to read_info
            read_info = add_entry_to_read_info(
                read_info=read_info,
                orig_read_id=read_id,
                new_read_id=read_id,
                readlen=len(read_seq),
                start=0,
                end=0,
                strand="=",
                adapter_config="ambiguous",
                lab="no_adapters",
            )
    return read_info


def write_stranded_fastq(tmp_fastq, read_info, single_end_mode):
    """Write stranded FASTQ entries based on read_info and return the output FASTQ file name."""
    tmp_stranded_fastq = tmp_fastq.replace(".fastq", ".stranded.fastq.gz")

    with pysam.FastxFile(tmp_fastq) as f_in:
        with gzip.open(tmp_stranded_fastq, "wb") as f_out:
            for entry in f_in:
                read_id = entry.name.split(" ")[0]
                if read_info.get(read_id):
                    for subread_id in read_info[read_id].keys():
                        d = read_info[read_id][subread_id]
                        if single_end_mode:  # TODO
                            subread_seq = str(entry.sequence)
                            subread_quals = entry.quality
                            if d["strand"] == "-":
                                subread_seq = subread_seq[::-1].translate(
                                    COMPLEMENT_TRANS
                                )
                                subread_quals = subread_quals[::-1]
                        else:
                            subread_seq = str(entry.sequence)
                            subread_quals = entry.quality

                            if d["strand"] == "-":
                                subread_seq = subread_seq[::-1].translate(
                                    COMPLEMENT_TRANS
                                )
                                subread_quals = subread_quals[::-1]
                        f_out.write(f"@{subread_id}\n".encode())
                        f_out.write(f"{subread_seq}\n".encode())
                        f_out.write(b"+\n")
                        f_out.write(f"{subread_quals}\n".encode())
                else:
                    pass

    return tmp_stranded_fastq


def count_reads(filename):
    """Count and return the number of reads in a FASTQ file, handling both gzipped and uncompressed files."""
    count = 0
    open_func = gzip.open if filename.endswith(".gz") else open
    with open_func(filename, "rb") as f:  # Open in binary mode
        buf = io.BufferedReader(f)
        for line in buf:
            count += 1
    return count // 4


def get_subread_info(read_info):
    """Extract and return subread information from the read_info dictionary as a list."""
    subread_info = []
    for read_id, subread_d in read_info.items():
        for subread_id, attr_d in subread_d.items():
            attr_d["read_id"] = subread_id
            subread_info.append(attr_d)
    return subread_info


def write_tmp_table(tmp_fastq, subread_info):
    """Write a temporary subread information table to disk and return the file name."""
    df = pd.DataFrame.from_records(subread_info)
    tmp_table = tmp_fastq.replace(".fastq.gz", ".info.tsv")
    df.to_csv(tmp_table, sep="\t", index=False)
    return tmp_table


def process_batch(tup):
    """
    Process a batch of reads from a temporary FASTQ file and return the processed FASTQ
    and table file names.

    Args:
        tup (tuple): A tuple containing the temporary FASTQ file name and the arguments object.

    Returns:
        tuple: A tuple containing the stranded temporary FASTQ file name and the temporary table file name.
    """
    tmp_fastq = tup[0]
    args = tup[1]

    # Read all entries from the temporary FASTQ file into a list
    with pysam.FastxFile(tmp_fastq) as f_in:
        reads = [entry for entry in f_in]

    # Parse the reads to find adapter alignments and determine valid adapter pairs
    read_info = parse_parasail(
        reads,
        single_end_mode=args.single_end_mode,
        adapter1_seq=args.adapter1_seq,
        adapter2_seq=args.adapter2_seq,
        min_adapter_match=args.min_adapter_match,
    )

    # Write the stranded FASTQ entries based on the parsed read information
    stranded_tmp_fastq = write_stranded_fastq(
        tmp_fastq, read_info, args.single_end_mode
    )

    # Extract subread information from the read_info dictionary
    subread_info = get_subread_info(read_info)

    # Write the subread information to a temporary table file
    tmp_table = write_tmp_table(stranded_tmp_fastq, subread_info)

    # Return the stranded FASTQ file name and the temporary table file name
    return stranded_tmp_fastq, tmp_table


def write_output_table(tmp_tables, args):
    """Merge temporary tables and write the final adapter configuration TSV output."""
    if len(tmp_tables) > 1:
        pd.concat([pd.read_csv(d, sep="\t") for d in tmp_tables], axis=0).to_csv(
            args.output_tsv, sep="\t", index=False
        )
    else:
        shutil.copy(tmp_tables[0], args.output_tsv)


def write_fq_out(tmp_fastqs, args):
    """Merge temporary FASTQ files into the final output file and remove the temporary directory."""
    with open(args.fq_out, "wb") as f_out:
        for tmp_fastq in tmp_fastqs:
            with open(tmp_fastq, "rb") as f_:
                shutil.copyfileobj(f_, f_out)

    shutil.rmtree(args.tempdir, ignore_errors=True)


def write_tmp_fastx_files_for_processing(n_batches, args):
    """Create temporary FASTQ files for each batch and return a dictionary of file names."""
    fastq_fns = {}
    fastq_fs = {}
    for i in range(1, n_batches + 1):
        tmp_prefix = f"tmp.chunk.{i}."
        tmp_fastq = tempfile.NamedTemporaryFile(
            prefix=tmp_prefix, suffix=".fastq", dir=args.tempdir, delete=False
        )
        fastq_fns[i] = tmp_fastq.name
        fastq_fs[i] = open(tmp_fastq.name, "w")

    batch_id = 1
    with pysam.FastxFile(args.fq_in, "r") as f_in:
        for i, entry in enumerate(f_in):
            if (i > 0) & (i % args.batch_size == 0):
                batch_id += 1
            f_out = fastq_fs[batch_id]
            f_out.write("@" + str(entry.name) + "\n")
            f_out.write(str(entry.sequence) + "\n")
            f_out.write(str("+\n"))
            f_out.write(str(entry.quality) + "\n")

    [f_out.close() for f_out in fastq_fs.values()]

    return fastq_fns


def main(args):
    """Main function to execute the processing pipeline, writing outputs, and reporting run time."""
    start_time = time.time()
    print(
        f"Run started at: {time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(start_time))}"
    )

    # Print run settings for log files
    print(
        f"adapter_scan version:             adapter_scan_parasail_v2\n"
        f"Input FASTQ:                      {args.fq_in}\n"
        f"Output FASTQ:                     {args.fq_out}\n"
        f"Output TSV:                       {args.output_tsv}\n"
        f"Threads:                          {args.threads}\n"
        f"Batch size:                       {args.batch_size}\n"
        f"Adapter1 sequence:                {args.adapter1_seq}\n"
        f"Adapter2 sequence:                {args.adapter2_seq}\n"
        f"Minimum adapter % match:          {args.min_adapter_match}\n"
    )

    # Check for identical adapter sequences or reverse complements
    if not args.single_end_mode:
        adapter1_set = set(args.adapter1_seq)
        adapter2_set = set(args.adapter2_seq)
        adapter2_rc_set = set(
            seq[::-1].translate(COMPLEMENT_TRANS) for seq in args.adapter2_seq
        )
        common_adapters = adapter1_set.intersection(adapter2_set).union(
            adapter1_set.intersection(adapter2_rc_set)
        )
        if common_adapters:
            print(
                f"Warning: The following adapter sequences are present in both adapter1 and adapter2 (including reverse complements): {common_adapters}"
            )
            print_log(f"Using single_end_mode\n")
            args.single_end_mode = True
        else:
            args.single_end_mode = False
    else:
        print("Single adapter mode enabled due to no adapter2 sequences provided\n")

    # Print all sequences used in the search
    print_log(f"All sequences used in the search:")
    for i, seq in enumerate(args.adapter1_seq):
        print(f"Adapter1-{i}: {seq}")
    print("")
    if not args.single_end_mode:
        for i, seq in enumerate(args.adapter2_seq):
            print(f"Adapter2-{i}: {seq}")
        print("")

    # If specified batch size is > total number of reads, reduce batch size
    print_log("Counting reads")
    n_reads = count_reads(args.fq_in)
    print_log(f"Total number of reads: {n_reads}")

    args.batch_size = min(n_reads, args.batch_size)
    n_batches = int(np.ceil(n_reads / args.batch_size))

    # Create temp directory
    if os.path.exists(args.tempdir):
        shutil.rmtree(args.tempdir, ignore_errors=True)
    os.mkdir(args.tempdir)

    fastq_fns = write_tmp_fastx_files_for_processing(n_batches, args)

    print_log(f"Processing {n_batches} batches of {args.batch_size} reads")
    func_args = [(fn, args) for fn in fastq_fns.values()]

    if args.threads == 1 or n_batches == 1:
        results = [process_batch(arg) for arg in func_args]
    else:
        with multiprocessing.Pool(args.threads) as p:
            results = list(p.imap(process_batch, func_args))

    tmp_fastqs, tmp_tables = zip(*results)

    # Merge temp tables and fastqs then clean up
    print_log(f"Writing output table to {args.output_tsv}")
    write_output_table(tmp_tables, args)

    print_log(f"Writing stranded fastq to {args.fq_out}")
    write_fq_out(tmp_fastqs, args)
    print("")

    end_time = time.time()
    print(
        f"Run finished at: {time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(end_time))}"
    )
    print(f"Total run time: {end_time - start_time:.2f} seconds")


if __name__ == "__main__":
    args = parse_args()
    main(args)
