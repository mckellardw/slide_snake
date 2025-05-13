import argparse
import pysam
import parasail
import time
import os
import sys

# Short-read usage:
## Visium
"""
python scripts/py/fastq_call_bc_umi_from_adapter_v2.py \
    --fq_in out/test_Vis_yPAP_3C/short_read/tmp/cut_R1.fq.gz \
    --tsv_out out/test_Vis_yPAP_3C/short_read/barcodes_umis/visium/read_barcodes.tsv \
    --bc_adapters 0 \
    --bc_lengths 16 \
    --bc_offsets 0 \
    --bc_positions right \
    --bc_min_adapter_match 2 \
    --umi_adapters 0 \
    --umi_lengths 12 \
    --umi_offsets 16 \
    --umi_positions right \
    --umi_min_adapter_match 2 \
    --threads 1
"""


# Long-read usage:
## SlideSeq
""" 
python scripts/py/fastq_call_bc_umi_from_adapter_v2.py \
    --fq_in data/test_seeker/ont/tmp/cut_R1.fq.gz \
    --tsv_out barcodes.tsv \
    --bc_adapters TCTTCAGCGTTCCCGAGA TCTTCAGCGTTCCCGAGA \
    --bc_lengths 8 6 \
    --bc_offsets 0 \
    --bc_positions left right \
    --bc_min_adapter_match 2 \
    --umi_adapters TCTTCAGCGTTCCCGAGA \
    --umi_lengths 7 \
    --umi_offsets 6 \
    --umi_positions right \
    --umi_min_adapter_match 2 \
    --threads 1
"""

## Visium
""" 
python scripts/py/fastq_call_bc_umi_from_adapter_v2.py \
    --fq_in data/test_visium/vy3C_1k_ont.fq.gz \
    --tsv_out sandbox/barcodes.tsv \
    --bc_adapters CTACACGACGCTCTTCCGATCT \
    --bc_lengths 16 \
    --bc_offsets 0 \
    --bc_positions right \
    --bc_min_adapter_match 2 \
    --umi_adapters CTACACGACGCTCTTCCGATCT \
    --umi_lengths 12 \
    --umi_offsets 16 \
    --umi_positions right \
    --umi_min_adapter_match 2 \
    --threads 1
"""

# microST - single-sector
"""
python scripts/py/fastq_call_bc_umi_from_adapter_v2.py \
    --fq_in out/E095_A_sub/tmp/ont/cut_R1.fq.gz \
    --tsv_out out/E095_A_sub/ont/barcodes_umis/microST_ssv1_matchLinker_total/read_barcodes.tsv \
    --bc_adapters TGATGCCACACTGA TGATGCCACACTGA \
    --bc_lengths 10 10 \
    --bc_offsets 0 0 \
    --bc_positions left right \
    --bc_min_adapter_match 2 \
    --umi_adapters CTACACGACGCTCTTCCGATCT \
    --umi_lengths 12 \
    --umi_offsets 0 \
    --umi_positions right \
    --umi_min_adapter_match 2 \
    --threads 1
"""


def currentTime():
    return time.strftime("%D | %H:%M:%S", time.localtime())


def print_log(msg):
    """Print a log message with the current time."""
    print(f"{currentTime()} | {msg}")


def print_error(msg):
    """Print an error message with the current time."""
    print(f"{currentTime()} | ERROR: {msg}", file=sys.stderr)
    sys.exit(1)


def parse_args():
    parser = argparse.ArgumentParser(
        description="Extract barcodes from FASTQ file and write them in the header."
    )
    parser.add_argument("--fq_in", required=True, help="Path to the input FASTQ file.")
    parser.add_argument(
        "--tsv_out", required=True, help="Path to the output FASTQ file."
    )
    parser.add_argument(
        "--bc_adapters",
        nargs="+",
        required=True,
        type=lambda x: int(x) if x.isdigit() else x,
        help="List of adapter sequences or positions used for barcode extraction.",
    )
    parser.add_argument(
        "--bc_positions",
        nargs="+",
        required=True,
        help="List of barcode positions (left or right).",
    )
    parser.add_argument(
        "--bc_lengths",
        nargs="+",
        type=int,
        required=True,
        help="List of barcode lengths.",
    )
    parser.add_argument(
        "--bc_offsets",
        nargs="+",
        type=int,
        required=True,
        help="List of offsets from adapter sequence. Keep in mind the `bc_position`",
    )
    parser.add_argument(
        "--bc_min_adapter_match",
        nargs="+",
        type=float,
        required=True,
        help="Fraction anchor sequence match required for BC calling. [0.7]",
    )
    parser.add_argument(
        "--umi_adapters",
        nargs="+",
        required=True,
        type=lambda x: int(x) if x.isdigit() else x,
        help="List of adapter sequences or positions used for UMI extraction.",
    )
    parser.add_argument(
        "--umi_positions",
        nargs="+",
        required=True,
        help="List of UMI positions (`left` or `right`).",
    )
    parser.add_argument(
        "--umi_lengths",
        nargs="+",
        type=int,
        required=True,
        help="List of UMI lengths.",
    )
    parser.add_argument(
        "--umi_offsets",
        nargs="+",
        type=int,
        required=True,
        help="List of offsets from adapter sequence. Keep in mind the `umi_position`",
    )
    parser.add_argument(
        "--umi_min_adapter_match",
        nargs="+",
        type=float,
        required=True,
        help="Fraction anchor sequence match required for UMI calling. [0.7]",
    )
    parser.add_argument(
        "--threads",
        type=int,
        required=True,
        help="Number of threads to use.",
    )
    parser.add_argument(
        "--stats_out",
        required=False,
        default=None,
        help="Path to the output TSV file for statistics.",
    )
    args = parser.parse_args()

    return args


def align_parasail(read, adapter, min_adapter_match=0.7, matrix=None):
    """
    Align sequences using parasail (Smith-Waterman local alignment)
    - source: https://github.com/jeffdaily/parasail-python
    """
    # Create scoring matrix and perform alignment
    if matrix is None:
        matrix = parasail.matrix_create("ACGT", 1, -1)

    # s1=query, s2=target
    alignment = parasail.sw_trace_striped_32(
        s1=adapter, s2=read, open=10, extend=1, matrix=matrix
    )

    # Check if alignment meets the minimum score threshold based on min_adapter_match
    min_score = round(min_adapter_match * len(adapter))
    if alignment.score >= min_score:
        start = alignment.end_ref - len(adapter) + 1
        end = alignment.end_ref + 1
        return alignment.score, start, end

    return None, None, None


def align_parasail_or_hardcoded(read, adapter, min_adapter_match=0.7, verbose=False):
    """
    Align sequences using parasail (Smith-Waterman local alignment) or,
    if an integer is passed for the adapter, use that as a hardcoded position.
    Supports both barcode and UMI extraction.
    """
    if isinstance(adapter, int):
        # Use hardcoded position for barcode or UMI extraction:
        start = adapter
        end = adapter
        return 420, start, end
    else:
        # Use parasail alignment
        return align_parasail(read, adapter, min_adapter_match)


# Simple python version of R rep() function
def rep(val, n):
    return [val for i in range(0, n)]


def main(
    fq_in,
    tsv_out,
    bc_adapters,
    bc_positions,
    bc_lengths,
    bc_offsets,
    bc_min_adapter_match,
    umi_adapters,
    umi_positions,
    umi_lengths,
    umi_offsets,
    umi_min_adapter_match,
):
    # Ensure the output directory exists
    os.makedirs(os.path.dirname(tsv_out), exist_ok=True)

    null_bc_string = "-"

    bc_match_count = 0
    bc_missing_count = 0
    no_bc_count = 0
    umi_match_count = 0
    umi_missing_count = 0
    no_umi_count = 0
    read_count = 0

    BC_RANGES = list(
        zip(bc_adapters, bc_positions, bc_lengths, bc_offsets, bc_min_adapter_match)
    )
    UMI_RANGES = list(
        zip(
            umi_adapters, umi_positions, umi_lengths, umi_offsets, umi_min_adapter_match
        )
    )

    with pysam.FastxFile(fq_in) as fastq:
        with open(tsv_out, "w") as outfile:
            for read in fastq:
                if read_count % 5000000 == 0:
                    print(f"{currentTime()} - {read_count} reads processed...")

                bc_to_write = []
                for adapter, position, length, offset, min_adapter_match in BC_RANGES:
                    if isinstance(adapter, int) or len(read.sequence) >= len(adapter):
                        align_score, start, end = align_parasail_or_hardcoded(
                            read=read.sequence,
                            adapter=adapter,
                            min_adapter_match=min_adapter_match,
                        )
                        # if align_score < 45: #debug
                        # print(alignment)
                        # print(align_score)
                        # print(read.sequence[start:end])

                        if align_score is not None:
                            # if align_score > 2 * len(adapter):
                            if position == "left":
                                bc_seq = read.sequence[
                                    (start - offset - length) : (start - offset)
                                ]
                            elif position == "right":
                                bc_seq = read.sequence[
                                    (end + offset) : (end + offset + length)
                                ]
                            else:
                                print(
                                    f"Incorrect barcode position [{position}] specified!"
                                )

                            # Don't write partial barcodes
                            if len(bc_seq) == length:
                                bc_to_write.append(bc_seq)
                            else:
                                bc_to_write.append(null_bc_string)
                        else:
                            bc_to_write.append(null_bc_string)
                # end BC alignment

                # Write barcode to file
                ## Only write BC if all positions are found
                if bc_to_write.count(null_bc_string) == 0:
                    bc_match_count += 1

                    # Proceed w/ UMI alignment
                    umi_to_write = []
                    for (
                        adapter,
                        position,
                        length,
                        offset,
                        min_adapter_match,
                    ) in UMI_RANGES:
                        if isinstance(adapter, int) or len(read.sequence) >= len(
                            adapter
                        ):
                            align_score, start, end = align_parasail_or_hardcoded(
                                read=read.sequence,
                                adapter=adapter,
                                min_adapter_match=min_adapter_match,
                            )

                            # if align_score > 2 * len(adapter):
                            if align_score is not None:
                                if position == "left":
                                    umi_seq = read.sequence[
                                        (start - offset - length) : (start - offset)
                                    ]
                                elif position == "right":
                                    umi_seq = read.sequence[
                                        (end + offset) : (end + offset + length)
                                    ]
                                else:
                                    print(
                                        f"Incorrect UMI position [{position}] specified!"
                                    )

                                # Don't write partial barcodes
                                if len(umi_seq) == length:
                                    umi_to_write.append(umi_seq)
                                else:
                                    umi_to_write.append(null_bc_string)
                            else:
                                umi_to_write.append(null_bc_string)
                    # end UMI alignment

                    outfile.write(
                        "{}\t{}\t{}\n".format(
                            read.name, "\t".join(bc_to_write), "\t".join(umi_to_write)
                        )
                    )
                    # Write UMI to file
                    if umi_to_write.count(null_bc_string) == 0:
                        umi_match_count += 1
                    elif umi_to_write.count(null_bc_string) == len(umi_adapters):
                        no_umi_count += 1
                    else:
                        umi_missing_count += 1
                elif bc_to_write.count(null_bc_string) == len(bc_adapters):
                    no_bc_count += 1
                else:
                    bc_missing_count += 1

                read_count += 1
            # end fq iterator
        # end tsv open
    # end fq open

    print("")
    print(f"{currentTime()} - {read_count:,} total reads processed.")

    return (
        read_count,
        bc_match_count,
        bc_missing_count,
        no_bc_count,
        umi_match_count,
        umi_missing_count,
        no_umi_count,
    )


if __name__ == "__main__":
    args = parse_args()

    # param checks ------
    if (
        len(args.bc_min_adapter_match) != len(args.bc_adapters)
        and len(args.bc_min_adapter_match) == 1
    ):
        args.bc_min_adapter_match = rep(
            val=args.bc_min_adapter_match[0], n=len(args.bc_adapters)
        )

    if len(args.bc_offsets) != len(args.bc_adapters) and len(args.bc_offsets) == 1:
        args.bc_offsets = rep(val=args.bc_offsets[0], n=len(args.bc_adapters))

    if (
        len(args.umi_min_adapter_match) != len(args.umi_adapters)
        and len(args.umi_min_adapter_match) == 1
    ):
        args.umi_min_adapter_match = rep(
            val=args.umi_min_adapter_match[0], n=len(args.umi_adapters)
        )

    if len(args.umi_offsets) != len(args.umi_adapters) and len(args.umi_offsets) == 1:
        args.umi_offsets = rep(val=args.umi_offsets[0], n=len(args.umi_adapters))

    # Set offsets to 0 for hardcoded adapter positions.
    for i, adapter in enumerate(args.bc_adapters):
        if isinstance(adapter, int):
            args.bc_offsets[i] = 0
    for i, adapter in enumerate(args.umi_adapters):
        if isinstance(adapter, int):
            args.umi_offsets[i] = 0

    # Print run settings for log files ----
    print(
        f"input fastq:                      {args.fq_in}\n"
        f"output tsv:                       {args.tsv_out}\n"
        f"BC adapter sequence(s):           {args.bc_adapters}\n"
        f"BC position(s):                   {args.bc_positions}\n"
        f"BC offset(s):                     {args.bc_offsets}\n"
        f"BC length(s):                     {args.bc_lengths}\n"
        f"BC adapter mismatches allowed:    {args.bc_min_adapter_match}\n"
        f"UMI adapter sequence(s):          {args.umi_adapters}\n"
        f"UMI position(s):                  {args.umi_positions}\n"
        f"UMI offset(s):                    {args.umi_offsets}\n"
        f"UMI length(s):                    {args.umi_lengths}\n"
        f"UMI adapter mismatches allowed:   {args.umi_min_adapter_match}\n"
        f"threads:                          {args.threads}\n"
    )

    # Print if hardcoded positions are being used
    if any(isinstance(adapter, int) for adapter in args.bc_adapters):
        print("Using hardcoded positions for barcodes.")
    if any(isinstance(adapter, int) for adapter in args.umi_adapters):
        print("Using hardcoded positions for UMIs.")

    print_log("Starting run...")
    print("")

    # Run main function ------
    if args.threads > 1:
        # # Source: https://superfastpython.com/multiprocessing-pool-for-loop/
        # import multiprocessing

        # # Split .fq file into {n_core} chunks
        # os.system(
        #     # Custom script to split .fq w/ sed in parallel
        #     f"""
        #     python scripts/py/splitNfqs.py {args.fq_in} {args.threads} {args.threads}
        #     """
        # )

        # # Trim each chunked file, in parallel
        # chunked_fqs_in = [
        #     f"{args.fq_in.strip('.fq.gz')}_{str(n).zfill(3)}.fq"
        #     for n in list(range(1, args.threads + 1))
        # ]

        # # multiprocessing
        # items = [
        #     (
        #         f"{args.fq_in.strip('.fq.gz')}_{str(n).zfill(3)}.fq",
        #         args.tsv_out,
        #         args.bc_adapters,
        #         args.bc_positions,
        #         args.bc_lengths,
        #         args.bc_offsets,
        #         args.bc_min_adapter_match,
        #         args.umi_adapters,
        #         args.umi_positions,
        #         args.umi_lengths,
        #         args.umi_offsets,
        #         args.umi_min_adapter_match,
        #     )
        #     for n in list(range(1, args.threads + 1))
        # ]

        # with multiprocessing.Pool(args.threads) as pool:
        #     multi_out = pool.starmap(main, items)

        # if os.path.isfile(args.fq_in.replace(".fq.gz", "_R2.fq.gz")):
        #     os.system(
        #         f"""
        #         rm {' '.join(chunked_fqs_in)}
        #         """
        #     )

        # TODO - need to split main() into multiple functions, to simplify writing output (without adding huge MEM overhead of returning entire output...)

        print(f"Multithreading not yet supported, running on 1 thread...")
        (
            read_count,
            bc_match_count,
            bc_missing_count,
            no_bc_count,
            umi_match_count,
            umi_missing_count,
            no_umi_count,
        ) = main(
            args.fq_in,
            args.tsv_out,
            args.bc_adapters,
            args.bc_positions,
            args.bc_lengths,
            args.bc_offsets,
            args.bc_min_adapter_match,
            args.umi_adapters,
            args.umi_positions,
            args.umi_lengths,
            args.umi_offsets,
            args.umi_min_adapter_match,
        )
    elif args.threads == 1:
        (
            read_count,
            bc_match_count,
            bc_missing_count,
            no_bc_count,
            umi_match_count,
            umi_missing_count,
            no_umi_count,
        ) = main(
            args.fq_in,
            args.tsv_out,
            args.bc_adapters,
            args.bc_positions,
            args.bc_lengths,
            args.bc_offsets,
            args.bc_min_adapter_match,
            args.umi_adapters,
            args.umi_positions,
            args.umi_lengths,
            args.umi_offsets,
            args.umi_min_adapter_match,
        )

    print_log("Run finished.")
    print(
        f"Summary:\n"
        f"   # Total reads processed:                  {read_count:,}\n"
        f"   # Reads w/ all barcodes found:            {bc_match_count:,} ({bc_match_count/read_count*100:.2f}%)\n"
        f"   # Reads w/ one or more barcodes missing:  {bc_missing_count:,} ({bc_missing_count/read_count*100:.2f}%)\n"
        f"   # Reads w/ no barcodes:                   {no_bc_count:,} ({no_bc_count/read_count*100:.2f}%)\n"
        f"\n"
        f"   # Reads w/ all UMIs found:                {umi_match_count:,} ({umi_match_count/read_count*100:.2f}%)\n"
        f"   # Reads w/ at least one UMI missing:      {umi_missing_count:,} ({umi_missing_count/read_count*100:.2f}%)\n"
        f"   # Reads w/ no UMIs:                       {no_umi_count:,} ({no_umi_count/read_count*100:.2f}%)\n"
        f"\n"
    )

    if args.stats_out:
        with open(args.stats_out, "w") as stats_file:
            stats_file.write(
                "total_reads\tbc_match_count\tbc_missing_count\tno_bc_count\tumi_match_count\tumi_missing_count\tno_umi_count\n"
                f"{read_count}\t{bc_match_count}\t{bc_missing_count}\t{no_bc_count}\t{umi_match_count}\t{umi_missing_count}\t{no_umi_count}\n"
            )
        print_log(f"Summary stats written to [{args.stats_out}]")
