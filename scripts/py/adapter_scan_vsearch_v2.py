import argparse
import gzip
import multiprocessing
import os
import pathlib
import shutil
import subprocess
import sys
import tempfile
import time

import numpy as np
import pandas as pd
import pysam

# Example usage:
"""
python scripts/py/adapter_scan_vsearch_v2.py \
    --fq_in "sandbox/adapter_scan/maxima/merged.fq.gz" \
    --fq_out "sandbox/adapter_scan/maxima/merged_stranded.fq.gz" \
    --output_tsv "sandbox/adapter_scan/maxima/adapter_scan.tsv" \
    --threads 56 \
    --batch_size 100000 \
    --adapter1_seq "CTACACGACGCTCTTCCGATCT" \
    --adapter2_seq "ATGTACTCTGCGTTGATACCACTGCTT" \
    --min_adapter_id 0.7 
"""


def parse_args():
    # Create argument parser
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
        help="Forward primer sequences (Read1). For TXG, Curio, etc. use [CTACACGACGCTCTTCCGATCT]",
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
        default=["ATGTACTCTGCGTTGATACCACTGCTT"],
        type=str,
    )

    parser.add_argument(
        "-i",
        "--min_adapter_id",
        help="Minimum adapter alignment identity for VSEARCH [0.7]",
        type=float,
        default=0.7,
    )

    parser.add_argument(
        "--only_strand_full_length",
        help="Do not try to strand-orient reads where either \
                        just a single adapter was found [False]",
        action="store_true",
        default=False,
    )

    parser.add_argument(
        "-a",
        "--adapters_fasta",
        help="Filename for adapter query sequences \
                        [adapter_seqs.fasta]",
        type=str,
        default="adapter_seqs.fasta",
    )

    parser.add_argument(
        "--verbosity",
        help="logging level: <=2 logs info, <=3 logs warnings",
        type=int,
        default=3,
    )

    # Parse arguments
    args = parser.parse_args()

    # Create temp dir and add that to the args object
    p = pathlib.Path(args.output_tsv)
    tempdir = tempfile.TemporaryDirectory(prefix="tmp.", dir=p.parents[0])
    args.tempdir = tempdir.name

    # Add tempdir path to the adapters_fasta
    args.adapters_fasta = os.path.join(args.tempdir, args.adapters_fasta)

    return args


# complement translation table with support for regex punctuation
COMPLEMENT_TRANS = str.maketrans(
    "ACGTWSMKRYBDHVNacgtwsmkrybdhvn", "TGCAWSKMYRVHDBNtgcawskmyrvhdbn"
)


def run_subprocess(cmd):
    """
    Run OS command and return stdout & stderr
    """
    p = subprocess.Popen(cmd, shell=True, stderr=subprocess.PIPE)
    stdout, stderr = p.communicate()
    return str(stdout), str(stderr)


def check_vsearch():
    stdout, stderr = run_subprocess("vsearch --quiet -h")
    if stderr.find("vsearch: command not found") > -1:
        print("Could not load find VSEARCH -- check installation")
        sys.exit(1)


def write_adapters_fasta(args):
    adapters = []
    for i, seq in enumerate(args.adapter1_seq):
        adapters.append(f">adapter1_f_{i}\n{seq}\n")
        adapters.append(f">adapter1_r_{i}\n{seq[::-1].translate(COMPLEMENT_TRANS)}\n")
    if not args.single_adapter_mode:
        for i, seq in enumerate(args.adapter2_seq):
            adapters.append(f">adapter2_f_{i}\n{seq}\n")
            adapters.append(
                f">adapter2_r_{i}\n{seq[::-1].translate(COMPLEMENT_TRANS)}\n"
            )
    with open(args.adapters_fasta, "w") as f:
        f.writelines(adapters)


def write_tmp_fasta(batch_reads, args):
    tmp_fasta = tempfile.NamedTemporaryFile(
        prefix="tmp.reads.", suffix=".fasta", dir=args.tempdir, delete=False
    )

    with open(tmp_fasta.name, mode="w") as f_out:
        for r in batch_reads:
            f_out.write(">" + str(r.name) + "\n")
            f_out.write(str(r.sequence) + "\n")
    return tmp_fasta.name


def call_vsearch(tmp_fastq, args):
    """ """
    # Convert batch FASTQ to FASTA for input into VSEARCH
    tmp_fasta = tmp_fastq.replace(".fastq", ".fasta")

    # Minimum adapter sequence length - adapters silently thrown out if shorter
    minseqlength = 15

    with pysam.FastxFile(tmp_fastq) as f_in:
        with open(tmp_fasta, "w") as f_out:
            for entry in f_in:
                f_out.write(">" + str(entry.name) + "\n")
                f_out.write(str(entry.sequence) + "\n")

    tmp_vsearch = tmp_fastq.replace(".fastq", ".vsearch.tsv")

    vsearch_cmd = "vsearch --usearch_global {fasta} --db {adapters} \
    --threads 1 --minseqlength {minseqlength} --maxaccepts 5 --id {id} --strand plus \
    --wordlength 3 --minwordmatches 10 --output_no_hits --userfields \
    'query+target+id+alnlen+mism+opens+qilo+qihi+qstrand+tilo+tihi+ql+tl' \
    --userout {output}".format(
        fasta=tmp_fasta,
        adapters=args.adapters_fasta,
        minseqlength=minseqlength,
        id=args.min_adapter_id,
        output=tmp_vsearch,
    )
    stdout, stderr = run_subprocess(vsearch_cmd)

    # Check if the output file exists
    if not os.path.exists(tmp_vsearch):
        raise FileNotFoundError(
            f"VSEARCH output file not found: {tmp_vsearch}. "
            f"Command: {vsearch_cmd}\nStdout: {stdout}\nStderr: {stderr}"
        )

    os.remove(tmp_fasta)
    return tmp_vsearch


def get_valid_adapter_pair_positions_in_read(read, args):
    valid_pairs_n = 0
    fl_pairs = []

    def write_valid_pair_dict(
        read, adapter_1_idx, adapter_2_idx, valid_pairs_n, strand
    ):
        read_id = read["query"].iloc[0]
        pair_str = (
            f"{read.iloc[adapter_1_idx]['target']}-{read.iloc[adapter_2_idx]['target']}"
        )
        fl_pair = {
            "read_id": "{}_{}".format(read_id, valid_pairs_n),
            "config": pair_str,
            "start": read.iloc[adapter_1_idx]["qilo"],
            "end": read.iloc[adapter_2_idx]["qihi"],
            "strand": strand,
        }
        return fl_pair, valid_pairs_n

    if args.single_adapter_mode:
        for i in range(len(args.adapter1_seq)):
            adapter1_f = f"adapter1_f_{i}"
            adapter1_r = f"adapter1_r_{i}"
            adapter_1_f_idxs = read.index[read["target"] == adapter1_f]
            adapter_1_r_idxs = read.index[read["target"] == adapter1_r]
            for adapter_1_f_idx in adapter_1_f_idxs:
                for adapter_1_r_idx in adapter_1_r_idxs:
                    if adapter_1_f_idx < adapter_1_r_idx:
                        strand = "+"
                        fl_pair, valid_pairs_n = write_valid_pair_dict(
                            read,
                            adapter_1_f_idx,
                            adapter_1_r_idx,
                            valid_pairs_n,
                            strand,
                        )
                        valid_pairs_n += 1
                        fl_pairs.append(fl_pair)
                    elif adapter_1_f_idx > adapter_1_r_idx:
                        strand = "-"
                        fl_pair, valid_pairs_n = write_valid_pair_dict(
                            read,
                            adapter_1_r_idx,
                            adapter_1_f_idx,
                            valid_pairs_n,
                            strand,
                        )
                        valid_pairs_n += 1
                        fl_pairs.append(fl_pair)
    else:
        compat_adapters = {
            f"adapter1_f_{i}": f"adapter2_f_{j}"
            for i in range(len(args.adapter1_seq))
            for j in range(len(args.adapter2_seq))
        }
        compat_adapters.update(
            {
                f"adapter1_r_{i}": f"adapter2_r_{j}"
                for i in range(len(args.adapter1_seq))
                for j in range(len(args.adapter2_seq))
            }
        )

        for adapter1 in compat_adapters.keys():
            adapter_1_idxs = read.index[read["target"] == adapter1]
            for adapter_1_idx in adapter_1_idxs:
                adapter_2_idx = adapter_1_idx + 1
                if adapter_2_idx < read.shape[0]:
                    if read.iloc[adapter_2_idx]["target"] == compat_adapters[adapter1]:
                        strand = "+" if "f" in adapter1 else "-"
                        fl_pair, valid_pairs_n = write_valid_pair_dict(
                            read, adapter_1_idx, adapter_2_idx, valid_pairs_n, strand
                        )
                        valid_pairs_n += 1
                        fl_pairs.append(fl_pair)
    return fl_pairs


def add_entry_to_read_info(
    read_info,
    orig_read_id,
    new_read_id,
    readlen,
    start,
    end,
    fl,
    stranded,
    strand,
    orig_adapter_config,
    adapter_config,
    lab,
):
    read_info[orig_read_id][new_read_id] = {
        "readlen": readlen,
        "start": start,
        "end": end,
        "fl": fl,
        "stranded": stranded,
        "orig_strand": strand,
        "orig_adapter_config": orig_adapter_config,
        "adapter_config": adapter_config,
        "lab": lab,
    }
    return read_info


def parse_vsearch(tmp_vsearch, args):
    colnames = [
        "query",
        "target",
        "id",
        "alnlen",
        "mism",
        "opens",
        "qilo",
        "qihi",
        "qstrand",
        "tilo",
        "tihi",
        "ql",
        "tl",
    ]

    df = pd.read_csv(tmp_vsearch, sep="\t", header=None, names=colnames)

    read_info = {}
    for read_id, read in df.groupby("query"):
        read = read.sort_values("qilo").reset_index()
        orig_adapter_config = "-".join(read["target"])
        orig_read_id = read["query"].iloc[0]
        read_info[orig_read_id] = {}

        fl_pairs = get_valid_adapter_pair_positions_in_read(read, args)
        if len(fl_pairs) > 0:
            strands = [fl_pair["strand"] for fl_pair in fl_pairs]
            if all(s == "+" for s in strands):
                strand = "+"
            elif all(s == "-" for s in strands):
                strand = "-"
            else:
                strand = "*"
            for fl_pair in fl_pairs:
                fl = True
                stranded = True
                read_id = fl_pair["read_id"]
                readlen = fl_pair["end"] - fl_pair["start"]
                start = fl_pair["start"]
                end = fl_pair["end"]
                adapter_config = fl_pair["config"]
                lab = "full_len"
                read_info = add_entry_to_read_info(
                    read_info,
                    orig_read_id,
                    read_id,
                    readlen,
                    start,
                    end,
                    fl,
                    stranded,
                    strand,
                    orig_adapter_config,
                    adapter_config,
                    lab,
                )
        else:
            read_id = "{}_0".format(read["query"].iloc[0])
            fl = False
            stranded = False
            strand = "*"
            readlen = read["ql"].iloc[0]
            start = 0
            end = readlen - 1
            adapter_config = "-".join(list(read["target"].values))
            if adapter_config in ["adapter2_r-adapter2_f", "adapter2_f-adapter2_r"]:
                lab = "double_adapter2"
            elif adapter_config in ["adapter1_r-adapter1_f", "adapter1_f-adapter1_r"]:
                lab = "double_adapter1"
            elif adapter_config in ["adapter2_f", "adapter2_r"]:
                lab = "single_adapter2"
                if not args.only_strand_full_length:
                    stranded = True
                    if adapter_config == "adapter2_f":
                        strand = "+"
                        start = 0
                        end = read.iloc[0]["qihi"]
                        readlen = end - start
                    elif adapter_config == "adapter2_r":
                        strand = "-"
                        start = read.iloc[0]["qilo"]
                        end = read.iloc[0]["ql"] - 1
                        readlen = end - start
            elif adapter_config in ["adapter1_f", "adapter1_r"]:
                lab = "single_adapter1"
                if not args.only_strand_full_length:
                    stranded = True
                    if adapter_config == "adapter1_f":
                        strand = "+"
                        start = read.iloc[0]["qilo"]
                        end = read.iloc[0]["ql"] - 1
                        readlen = end - start
                    elif adapter_config == "adapter1_r":
                        strand = "-"
                        start = 0
                        end = read.iloc[0]["qihi"]
                        readlen = end - start
            elif adapter_config == "*":
                lab = "no_adapters"
            else:
                lab = "other"
            read_info = add_entry_to_read_info(
                read_info,
                orig_read_id,
                read_id,
                readlen,
                start,
                end,
                fl,
                stranded,
                strand,
                orig_adapter_config,
                adapter_config,
                lab,
            )
    return read_info, colnames


def revcomp_adapter_config(adapters_string):
    d = {
        "adapter1_f": "adapter1_r",
        "adapter1_r": "adapter1_f",
        "adapter2_f": "adapter2_r",
        "adapter2_r": "adapter2_f",
    }

    # Handle the new adapter identifiers
    def get_revcomp_adapter(adapter):
        parts = adapter.split("_")
        if len(parts) == 3:
            return f"{parts[0]}_{d[parts[0] + '_' + parts[1]]}_{parts[2]}"
        return adapter

    rc_string = "-".join(
        [get_revcomp_adapter(a) for a in adapters_string.split("-")[::-1]]
    )
    return rc_string


def write_stranded_fastq(tmp_fastq, read_info, args):
    tmp_stranded_fastq = tmp_fastq.replace(".fastq", ".stranded.fastq.gz")

    with pysam.FastxFile(tmp_fastq) as f_in:
        with gzip.open(tmp_stranded_fastq, "wb") as f_out:
            for entry in f_in:
                read_id = entry.name.split(" ")[0]
                if read_info.get(read_id):
                    for subread_id in read_info[read_id].keys():
                        d = read_info[read_id][subread_id]
                        subread_seq = str(entry.sequence[d["start"] : d["end"]])
                        subread_quals = entry.quality[d["start"] : d["end"]]
                        if d["orig_strand"] == "-":
                            rc_config = revcomp_adapter_config(d["adapter_config"])
                            d["adapter_config"] = rc_config
                            subread_seq = subread_seq[::-1].translate(COMPLEMENT_TRANS)
                            subread_quals = subread_quals[::-1]
                        f_out.write(f"@{subread_id}\n".encode())
                        f_out.write(f"{subread_seq}\n".encode())
                        f_out.write(b"+\n")
                        f_out.write(f"{subread_quals}\n".encode())
                else:
                    pass

    return tmp_stranded_fastq


def open_fastq(fastq):
    if fastq.split(".")[-1] == "gz":
        f = gzip.open(fastq, "rt")
    else:
        f = open(fastq)
    return f


def count_reads(fastq):
    number_lines = 0
    with open_fastq(fastq) as f:
        for line in f:
            number_lines += 1
    f.close()
    return number_lines / 4


def batch_iterator(iterator, args):
    """Returns lists of length batch_size.

    This can be used on any iterator, for example to batch up
    SeqRecord objects from Bio.SeqIO.parse(...), or to batch
    Alignment objects from Bio.AlignIO.parse(...), or simply
    lines from a file handle.

    This is a generator function, and it returns lists of the
    entries from the supplied iterator.  Each list will have
    batch_size entries, although the final list may be shorter.
    """
    entry = True  # Make sure we loop once
    while entry:
        batch = []
        while len(batch) < args.batch_size:
            try:
                entry = next(iterator)
            except StopIteration:
                entry = None
            if entry is None:
                # End of file
                break
            batch.append(entry)
        if batch:
            yield batch, args


def get_subread_info(read_info):
    subread_info = []
    for read_id, subread_d in read_info.items():
        for subread_id, attr_d in subread_d.items():
            attr_d["read_id"] = subread_id
            subread_info.append(attr_d)
    return subread_info


def write_tmp_table(tmp_fastq, subread_info):
    df = pd.DataFrame.from_records(subread_info)
    tmp_table = tmp_fastq.replace(".fastq.gz", ".info.tsv")
    df.to_csv(tmp_table, sep="\t", index=False)
    return tmp_table


def process_batch(tup):
    tmp_fastq = tup[0]
    args = tup[1]

    tmp_vsearch = call_vsearch(tmp_fastq, args)
    read_info, vsearch_cols = parse_vsearch(tmp_vsearch, args)
    stranded_tmp_fastq = write_stranded_fastq(tmp_fastq, read_info, args)
    subread_info = get_subread_info(read_info)
    tmp_table = write_tmp_table(stranded_tmp_fastq, subread_info)

    return stranded_tmp_fastq, tmp_table, vsearch_cols


def write_output_table(tmp_tables, args):
    """ """
    if len(tmp_tables) > 1:
        pd.concat([pd.read_csv(d, sep="\t") for d in tmp_tables], axis=0).to_csv(
            args.output_tsv, sep="\t", index=False
        )
    else:
        shutil.copy(tmp_tables[0], args.output_tsv)


def write_fq_out(tmp_fastqs, args):
    """ """
    with open(args.fq_out, "wb") as f_out:
        for tmp_fastq in tmp_fastqs:
            with open(tmp_fastq, "rb") as f_:
                shutil.copyfileobj(f_, f_out)

    shutil.rmtree(args.tempdir, ignore_errors=True)


def write_tmp_fastx_files_for_processing(n_batches, args):
    write_adapters_fasta(args)

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
    start_time = time.time()
    print(
        f"Run started at: {time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(start_time))}"
    )

    check_vsearch()

    # Print run settings for log files
    print(
        f"Input FASTQ:                      {args.fq_in}\n"
        f"Output FASTQ:                     {args.fq_out}\n"
        f"Output TSV:                       {args.output_tsv}\n"
        f"Threads:                          {args.threads}\n"
        f"Batch size:                       {args.batch_size}\n"
        f"Adapter1 sequence:                {args.adapter1_seq}\n"
        f"Adapter2 sequence:                {args.adapter2_seq}\n"
        f"Minimum adapter identity:         {args.min_adapter_id}\n"
        f"Only strand full length:          {args.only_strand_full_length}\n"
        f"Adapters FASTA:                   {args.adapters_fasta}\n"
        f"Verbosity:                        {args.verbosity}\n"
    )

    # Check for identical adapter sequences or reverse complements
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
        print(f"Using single_adapter_mode")
        args.single_adapter_mode = True
    else:
        args.single_adapter_mode = False

    # Print all sequences used in the search
    print(f"All sequences used in the search:")
    for seq in args.adapter1_seq:
        print(f"Adapter1: {seq}")
    if not args.single_adapter_mode:
        for seq in args.adapter2_seq:
            print(f"Adapter2: {seq}")

    # If specified batch size is > total number of reads, reduce batch size
    print("Counting reads")
    n_reads = count_reads(args.fq_in)
    args.batch_size = min(n_reads, args.batch_size)
    n_batches = int(np.ceil(n_reads / args.batch_size))

    # Create temp directory
    if os.path.exists(args.tempdir):
        shutil.rmtree(args.tempdir, ignore_errors=True)
    os.mkdir(args.tempdir)

    fastq_fns = write_tmp_fastx_files_for_processing(n_batches, args)

    print(f"Processing {n_batches} batches of {args.batch_size} reads")
    func_args = [(fn, args) for fn in fastq_fns.values()]

    if args.threads == 1 or n_batches == 1:
        results = [process_batch(arg) for arg in func_args]
    else:
        with multiprocessing.Pool(args.threads) as p:
            results = list(p.imap(process_batch, func_args))

    tmp_fastqs, tmp_tables, vsearch_cols = zip(*results)
    vsearch_cols = vsearch_cols[0]

    # Merge temp tables and fastqs then clean up
    print(f"Writing output table to {args.output_tsv}")
    write_output_table(tmp_tables, args)

    print(f"Writing stranded fastq to {args.fq_out}")
    write_fq_out(tmp_fastqs, args)

    end_time = time.time()
    print(
        f"Run finished at: {time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(end_time))}"
    )
    print(f"Total run time: {end_time - start_time:.2f} seconds")


if __name__ == "__main__":
    args = parse_args()
    main(args)
