import argparse
import gzip
import multiprocessing
import os
import pathlib
import shutil
import sys
import tempfile
import time

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
    --min_adapter_id 0.7 
"""



def currentTime():
    return time.strftime("%D | %H:%M:%S", time.localtime())

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
        nargs='+',
        default=["CTACACGACGCTCTTCCGATCT"],
        type=str,
    )

    parser.add_argument(
        "-r",
        "--adapter2_seq",
        help="Reverse complement of the primer sequences (rcTSO, minus the CCC prefix). \
            For TXG/Curio chemistries, use [ATGTACTCTGCGTTGATACCACTGCTT] (default). \
            For uMRT chemistries, use [GAGAGAGGAAAGAGAGAGAGAGGG]",
        nargs='+',
        default=["ATGTACTCTGCGTTGATACCACTGCTT"],
        type=str,
    )

    parser.add_argument(
        "-i",
        "--min_adapter_id",
        help="Minimum adapter alignment identity for parasail [0.7]",
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

    return args

# complement translation table with support for regex punctuation
COMPLEMENT_TRANS = str.maketrans(
    "ACGTWSMKRYBDHVNacgtwsmkrybdhvn", "TGCAWSKMYRVHDBNtgcawskmyrvhdbn"
)

def align_parasail(read, adapter, min_adapter_id):
    """
    Align sequences using parasail (Smith-Waterman local alignment)
    - source: https://github.com/jeffdaily/parasail-python
    """
    matrix = parasail.matrix_create("ACGT", 1, -1)
    result = parasail.sw_trace_striped_32(
        s1=adapter, s2=read, open=10, extend=1, matrix=matrix
    )

    mismatches = round((1 - min_adapter_id) * len(adapter))
    if result.score >= (len(adapter) - mismatches):
        return result
    return None


def get_valid_adapter_pair_positions_in_read(read, args):
    valid_pairs_n = 0
    fl_pairs = []

    def write_valid_pair_dict(read, adapter_1_idx, adapter_2_idx, valid_pairs_n, strand):
        read_id = read["query"].iloc[0]
        pair_str = (
            f"{read.iloc[adapter_1_idx]['target']}-{read.iloc[adapter_2_idx]['target']}"
        )
        fl_pair = {
            "read_id": "{}_{}".format(read_id, valid_pairs_n),
            "config": pair_str,
            "start": read.iloc[adapter_1_idx]["qilo"],
            "end": read.iloc[adapter_2_idx]["qihi"],
            "strand": strand
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
                            read, adapter_1_f_idx, adapter_1_r_idx, valid_pairs_n, strand
                        )
                        valid_pairs_n += 1
                        fl_pairs.append(fl_pair)
                    elif adapter_1_f_idx > adapter_1_r_idx:
                        strand = "-"
                        fl_pair, valid_pairs_n = write_valid_pair_dict(
                            read, adapter_1_r_idx, adapter_1_f_idx, valid_pairs_n, strand
                        )
                        valid_pairs_n += 1
                        fl_pairs.append(fl_pair)
    else:
        compat_adapters = {f"adapter1_f_{i}": f"adapter2_f_{j}" for i in range(len(args.adapter1_seq)) for j in range(len(args.adapter2_seq))}
        compat_adapters.update({f"adapter1_r_{i}": f"adapter2_r_{j}" for i in range(len(args.adapter1_seq)) for j in range(len(args.adapter2_seq))})

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

def parse_parasail(reads, args):
    read_info = {}
    for read in reads:
        read_id = read.name
        read_seq = read.sequence
        read_info[read_id] = {}
        alignments = []
        for i, adapter in enumerate(args.adapter1_seq):
            result = align_parasail(read_seq, adapter, args.min_adapter_id)
            if result:
                alignments.append({
                    "query": read_id,
                    "target": f"adapter1_f_{i}",
                    "qilo": result.end_query - len(adapter) + 1,
                    "qihi": result.end_query + 1,
                    "qstrand": "+",
                })
                result_rc = align_parasail(read_seq, adapter[::-1].translate(COMPLEMENT_TRANS), args.min_adapter_id)
                if result_rc:
                    alignments.append({
                        "query": read_id,
                        "target": f"adapter1_r_{i}",
                        "qilo": result_rc.end_query - len(adapter) + 1,
                        "qihi": result_rc.end_query + 1,
                        "qstrand": "-",
                    })
        if not args.single_adapter_mode:
            for i, adapter in enumerate(args.adapter2_seq):
                result = align_parasail(read_seq, adapter, args.min_adapter_id)
                if result:
                    alignments.append({
                        "query": read_id,
                        "target": f"adapter2_f_{i}",
                        "qilo": result.end_query - len(adapter) + 1,
                        "qihi": result.end_query + 1,
                        "qstrand": "+",
                    })
                    result_rc = align_parasail(read_seq, adapter[::-1].translate(COMPLEMENT_TRANS), args.min_adapter_id)
                    if result_rc:
                        alignments.append({
                            "query": read_id,
                            "target": f"adapter2_r_{i}",
                            "qilo": result_rc.end_query - len(adapter) + 1,
                            "qihi": result_rc.end_query + 1,
                            "qstrand": "-",
                        })
        if alignments:
            df = pd.DataFrame(alignments)
            df = df.sort_values("qilo").reset_index()
            orig_adapter_config = "-".join(df["target"])
            fl_pairs = get_valid_adapter_pair_positions_in_read(df, args)
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
                    new_read_id = fl_pair["read_id"]
                    readlen = fl_pair["end"] - fl_pair["start"]
                    start = fl_pair["start"]
                    end = fl_pair["end"]
                    adapter_config = fl_pair["config"]
                    lab = "full_len"
                    read_info = add_entry_to_read_info(
                        read_info,
                        read_id,
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
                    )
            else:
                new_read_id = f"{read_id}_0"
                fl = False
                stranded = False
                strand = "*"
                readlen = len(read_seq)
                start = 0
                end = readlen - 1
                adapter_config = "-".join(list(df["target"].values))
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
                            end = df.iloc[0]["qihi"]
                            readlen = end - start
                        elif adapter_config == "adapter2_r":
                            strand = "-"
                            start = df.iloc[0]["qilo"]
                            end = len(read_seq) - 1
                            readlen = end - start
                elif adapter_config in ["adapter1_f", "adapter1_r"]:
                    lab = "single_adapter1"
                    if not args.only_strand_full_length:
                        stranded = True
                        if adapter_config == "adapter1_f":
                            strand = "+"
                            start = df.iloc[0]["qilo"]
                            end = len(read_seq) - 1
                            readlen = end - start
                        elif adapter_config == "adapter1_r":
                            strand = "-"
                            start = 0
                            end = df.iloc[0]["qihi"]
                            readlen = end - start
                elif adapter_config == "*":
                    lab = "no_adapters"
                else:
                    lab = "other"
                read_info = add_entry_to_read_info(
                    read_info,
                    read_id,
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
                )
    return read_info

def revcomp_adapter_config(adapters_string):
    d = {
        "adapter1_f": "adapter1_r",
        "adapter1_r": "adapter1_f",
        "adapter2_f": "adapter2_r",
        "adapter2_r": "adapter2_f",
    }

    def get_revcomp_adapter(adapter):
        parts = adapter.split("_")
        if len(parts) == 3:
            return f"{parts[0]}_{d[parts[0] + '_' + parts[1]]}_{parts[2]}"
        return adapter

    rc_string = "-".join([get_revcomp_adapter(a) for a in adapters_string.split("-")[::-1]])
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
    entry = True
    while entry:
        batch = []
        while len(batch) < args.batch_size:
            try:
                entry = next(iterator)
            except StopIteration:
                entry = None
            if entry is None:
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

    with pysam.FastxFile(tmp_fastq) as f_in:
        reads = [entry for entry in f_in]
    read_info = parse_parasail(reads, args)
    stranded_tmp_fastq = write_stranded_fastq(tmp_fastq, read_info, args)
    subread_info = get_subread_info(read_info)
    tmp_table = write_tmp_table(stranded_tmp_fastq, subread_info)

    return stranded_tmp_fastq, tmp_table

def write_output_table(tmp_tables, args):
    if len(tmp_tables) > 1:
        pd.concat([pd.read_csv(d, sep="\t") for d in tmp_tables], axis=0).to_csv(
            args.output_tsv, sep="\t", index=False
        )
    else:
        shutil.copy(tmp_tables[0], args.output_tsv)

def write_fq_out(tmp_fastqs, args):
    with open(args.fq_out, "wb") as f_out:
        for tmp_fastq in tmp_fastqs:
            with open(tmp_fastq, "rb") as f_:
                shutil.copyfileobj(f_, f_out)

    shutil.rmtree(args.tempdir, ignore_errors=True)

def write_tmp_fastx_files_for_processing(n_batches, args):
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
    print(f"Run started at: {time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(start_time))}")

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
        f"Verbosity:                        {args.verbosity}\n"
    )

    # Check for identical adapter sequences or reverse complements
    adapter1_set = set(args.adapter1_seq)
    adapter2_set = set(args.adapter2_seq)
    adapter2_rc_set = set(seq[::-1].translate(COMPLEMENT_TRANS) for seq in args.adapter2_seq)
    common_adapters = adapter1_set.intersection(adapter2_set).union(adapter1_set.intersection(adapter2_rc_set))
    if common_adapters:
        print(f"Warning: The following adapter sequences are present in both adapter1 and adapter2 (including reverse complements): {common_adapters}")
        print(f"Using single_adapter_mode\n")
        args.single_adapter_mode = True
    else:
        args.single_adapter_mode = False

    # Print all sequences used in the search
    print(f"All sequences used in the search:")
    for seq in args.adapter1_seq:
        print(f"Adapter1: {seq}")
    print("")
    if not args.single_adapter_mode:
        for seq in args.adapter2_seq:
            print(f"Adapter2: {seq}")
        print("")

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

    tmp_fastqs, tmp_tables = zip(*results)

    # Merge temp tables and fastqs then clean up
    print(f"Writing output table to {args.output_tsv}")
    write_output_table(tmp_tables, args)

    print(f"Writing stranded fastq to {args.fq_out}")
    write_fq_out(tmp_fastqs, args)

    end_time = time.time()
    print(f"Run finished at: {time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(end_time))}")
    print(f"Total run time: {end_time - start_time:.2f} seconds")

if __name__ == "__main__":
    args = parse_args()
    main(args)
