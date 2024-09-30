# Script to remove the adapter sequence from R1

# TODO: add in read filtering for both R1 and R2 for reads which have alignment[0].score < 50

# Usage examples:
#   Snakemake:
#       python scripts/internal_adapter_trim_R1.py {params.INTERNAL_ADAPTER} {log} {threads} {params.TMPDIR} {input.MERGED_R1_FQ} {output.FINAL_R1_FQ}
#
#   bash:
#       python scripts/py/fastq_internal_adapter_trim_R1.py \
#           --adapter_seq TCTTCAGCGTTCCCGAGA \
#           --n_cores 8 \
#           --tmp_dir ./tmp/ \
#           --fq1_in data/test_seeker/heart_ctrl_1k_R1.fq.gz \
#           --fq1_out heart_ctrl_1k_R1_trimmed.fq.gz \
#           --min_adapter_start_pos 8

# Curio R1 structure:
#  [BB1- 8bp][Adapter Seq - 18bp][BB2|6bp][UMI-7bp][TTTTTTTTTTT]

# imports
import argparse
import sys
import os
import gzip
import time
import pysam
import parasail


def currentTime():
    return time.strftime("%D | %H:%M:%S", time.localtime())


def parse_args():
    parser = argparse.ArgumentParser(
        description="Script to remove the adapter sequence from R1."
    )
    parser.add_argument(
        "--adapter_seq", type=str, help="The internal adapter sequence."
    )
    parser.add_argument("--n_cores", type=int, help="Number of threads to use.")
    parser.add_argument("--tmp_dir", type=str, help="Directory for temporary files.")
    parser.add_argument("--fq1_in", type=str, help="Path to the merged R1 FASTQ file.")
    parser.add_argument(
        "--fq1_out", type=str, help="Path to the output final R1 FASTQ file."
    )
    parser.add_argument(
        "--min_adapter_start_pos",
        type=int,
        default=8,
        help="Minimum adapter start position.",
    )
    parser.add_argument(
        "--min_align_score", type=float, default=58, help="Minimum alignment score."
    )
    args = parser.parse_args()

    return args


def align_parasail(read, adapter, min_align_score=58, verbose=False):
    """
    Align sequences using parasail (Smith-Waterman local alignment)
    - source: https://github.com/jeffdaily/parasail-python
    """

    # Create a simple identity matrix (match = 1, mismatch = 0)
    matrix = parasail.matrix_create(alphabet="ACGT", match=4, mismatch=-0.5)
    alignment = parasail.sw(s1=adapter, s2=read, open=-6, extend=-6, matrix=matrix)

    # Check if alignment meets the minimum score threshold based on mismatches
    if alignment.score <= min_align_score:
        start = alignment.end_ref - len(adapter) + 1
        end = alignment.end_ref + 1
        return alignment.score, start, end
    return None, None, None


# Function to run internal trimming on a single .fq.gz file
def trim_fq(fq_in, fq_out, adapter_seq, min_adapter_start_pos, min_align_score):
    # Tallies for log file
    read_count = 0  # Read count
    ins_count = 0  # Insertion counter for BC_1
    del_count = 0  # Deletion counter for BC_1
    no_adapter_count = 0  # Tally of reads that are removed b/c missing adapter

    with gzip.open(fq_out, "wt") as out_handle:
        # Open the FASTQ file using pysam
        with pysam.FastxFile(fq_in) as fq:
            for read in fq:
                read_count += 1

                # Perform pairwise alignment to find the sequence with allowed score
                align_score, start, end = align_parasail(
                    read=read.sequence,
                    adapter=adapter_seq,
                    min_align_score=min_align_score,
                )

                # Account for reads with deletions in `BC_1`
                if start < min_adapter_start_pos:  # Deletion in BC_1
                    offset = min_adapter_start_pos - start
                    seq_out = (
                        "N" * offset
                        + read.sequence[start:min_adapter_start_pos]
                        + read.sequence[end:]
                    )
                    qual_out = (
                        "!" * offset
                        + read.quality[start:min_adapter_start_pos]
                        + read.quality[end:]
                    )

                    del_count += 1
                else:
                    if start > min_adapter_start_pos:  # Insertion in BC_1
                        ins_count += 1

                    ## Trim the base closest to adapter
                    seq_out = (
                        read.sequence[0:min_adapter_start_pos] + read.sequence[end:]
                    )
                    qual_out = (
                        read.quality[0:min_adapter_start_pos] + read.quality[end:]
                    )

                # Broken read; erase R1 and add `N` with qval=0 ('!')
                # Alignment score cutoff was determined by spot-checking ~100 bead barcodes, based on predicted adapter location
                # Seeker recommendation is min_align_score=58
                # TODO fix hardcode '22'
                if len(seq_out) < 22 or align_score < min_align_score:
                    seq_out = "N"
                    qual_out = "!"
                    no_adapter_count += 1

                # Write to new .fq.gz file
                out_handle.write(f"@{read.name}\n{seq_out}\n+\n{qual_out}\n")

    return [read_count, ins_count, del_count, no_adapter_count]


def main(args):
    if args.n_cores > 1:
        # Source: https://superfastpython.com/multiprocessing-pool-for-loop/
        import multiprocessing

        # Split .fq file into {n_core} chunks
        # TODO- move this code to be internal
        os.system(
            f""" 
            python scripts/py/splitNfqs.py {args.fq1_in} {args.n_cores} {args.n_cores}
            """
        )

        # Trim each chunked file, in parallel
        # tmp_fqs_out = [
        #     f"{args.fq1_in.replace('.fq.gz','')}_{str(n).zfill(3)}_trimmed.fq"
        #     for n in list(range(1, args.n_cores + 1))
        # ]
        items = [
            (
                f"{args.fq1_in.replace('.fq.gz','')}_{str(n).zfill(3)}.fq",
                f"{args.fq1_in.replace('.fq.gz','')}_{str(n).zfill(3)}_trimmed.fq",
                args.adapter_seq,
                args.min_adapter_start_pos,
                args.min_align_score,
            )
            for n in list(range(1, args.n_cores + 1))
        ]

        with multiprocessing.Pool(args.n_cores) as pool:
            multi_out = pool.starmap(trim_fq, items)

        read_count = 0
        ins_count = 0
        del_count = 0
        no_adapter_count = 0
        for i in list(range(0, len(multi_out))):
            read_count += multi_out[i][0]
            ins_count += multi_out[i][1]
            del_count += multi_out[i][2]
            no_adapter_count += multi_out[i][3]

        # Merge and compress chunked/trimmed fastqs
        os.system(
            f"""
            cat {args.fq1_in.replace('.fq.gz','')}_*_trimmed.fq > {args.fq1_out.replace('.gz','')}
            pigz -p{args.n_cores} {args.fq1_out.replace('.gz','')}
            """
        )

        # Remove tmp fastqs
        if os.path.isfile(args.fq1_out):
            os.system(
                f"""
                rm {args.fq1_in.replace('.fq.gz','')}_[0-9][0-9][0-9].fq {args.fq1_in.replace('.fq.gz','')}_*_trimmed.fq
                """
            )

        print(f"Total read count:         {read_count:,}")
        print(f"Insertion count in BC_1:  {ins_count:,}")
        print(f"Deletion count in BC_1:   {del_count:,}")
        print(f"Reads missing adapter: {no_adapter_count:,}")
    else:
        out = trim_fq(args.fq1_in, args.fq1_out, True)

        read_count, ins_count, del_count, no_adapter_count = out

        # Compress trimmed fastq
        os.system(
            f"""
            pigz -f -p{args.n_cores} {args.fq1_out.replace('.gz','')}
            """
        )

        print(f"Total read count:         {read_count:,}")
        print(f"Insertion count in BC_1:  {ins_count:,}")
        print(f"Deletion count in BC_1:   {del_count:,}")
        print(
            f"Reads trimmed below {args.min_adapter_start_pos}bp: {no_adapter_count:,}"
        )


if __name__ == "__main__":
    args = parse_args()

    # Run params
    print(f"Adapter sequence:                   {args.adapter_seq}")
    print(f"Input fastq:                        {args.fq1_in}")
    print(f"Output fastq:                       {args.fq1_out}")
    print(f"Minimum adapter start position:     {args.min_adapter_start_pos}")
    print(f"Minimum adapter alignment score:    {args.min_align_score}")
    print(f"{currentTime()} - Run starting...")

    main(args)
    print(f"{currentTime()} - Run complete.")
