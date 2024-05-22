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

# R1 structure:
#  [BB1- 8bp][Adapter Seq - 18bp][BB2|6bp][UMI-7bp][TTTTTTTTTTT]

# imports
import argparse
import sys
import os
import gzip
from Bio.SeqIO.QualityIO import FastqGeneralIterator
from Bio import Align


def init_aligner():
    # Useful resource: https://www.bioinformaticscrashcourse.com/10.1_Alignment.html
    # Curio data alignment parameters:
    ## match score = 4, mismatch = -0.5
    ## gap opening = -6, gap extension = -6
    aligner = Align.PairwiseAligner()
    aligner.mode = "local"  # Use 'local' for local alignment
    aligner.match_score = 4  # Match score
    aligner.mismatch_score = -0.5  # Mismatch score
    aligner.open_gap_score = -6  # Gap opening penalty
    aligner.extend_gap_score = -6  # Gap extension penalty

    return aligner


# Function to run internal trimming on a single .fq.gz file
def trim_fq(fq_in, fq_out, adapter_seq, min_adapter_start_pos, min_align_score):
    # Tallies for log file
    read_count = 0  # Read count
    ins_count = 0  # Insertion counter for BB_1
    del_count = 0  # Deletion counter for BB_1
    no_adapter_count = 0  # Tally of reads that are removed b/c missing adapter

    out_handle = open(fq_out, "w")

    if fq_in.endswith(".gz"):
        fq_iterator = FastqGeneralIterator(gzip.open(fq_in, "rt"))
    else:
        fq_iterator = FastqGeneralIterator(open(fq_in, "r"))

    # initialize aligner
    aligner = init_aligner()

    for title, seq, qual in fq_iterator:
        read_count += 1
        # Perform pairwise alignment to find the sequence with allowed mismatches
        alignments = aligner.align(seq, adapter_seq)

        # Get best alignment
        alignment = alignments[0] if alignments else None

        # For troubleshooting:
        # if (
        #     alignments[0].score < 60
        # ):  # and alignments[0].score > 55:# and alignments[0].start > 9:
        #     print(alignments[0])
        #     print(alignment.score)
        #     print(alignment.aligned[0][0][0])
        #     print(pairwise2.format_alignment(*alignments[0]))
        #     print(seq[0:8]+seq[alignments[0].end:])

        # Trim seq and qual
        start = alignment.aligned[0][0][0]
        end = alignment.aligned[0][0][1]

        # Acount for reads with deletions in `BB_1`
        if start < min_adapter_start_pos:  # Deletion in BB_1
            offset = min_adapter_start_pos - start
            seq_out = "N" * offset + seq[start:min_adapter_start_pos] + seq[end:]
            qual_out = "!" * offset + qual[start:min_adapter_start_pos] + qual[end:]

            del_count += 1
        else:
            if start > min_adapter_start_pos:  # Insertion in BB_1
                ins_count += 1

            ## Trim the base closest to adapter
            seq_out = seq[0:min_adapter_start_pos] + seq[end:]
            qual_out = qual[0:min_adapter_start_pos] + qual[end:]

        # Broken read; erase R1 and add `N` with qval=0 ('!')
        # Alignment score cuttoff was determined by spot-checking ~100 bead barcodes, based on predicted adapter location
        # Seeker recommendation is min_align_score=58
        # TODO fix hardcode '22'
        if len(seq_out) < 22 or alignment.score < min_align_score:
            seq_out = "N"
            qual_out = "!"
            no_adapter_count += 1

        # Write to new .fq.gz file
        out_handle.write(f"@{title}\n{seq_out}\n+\n{qual_out}\n")

    out_handle.close()

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
        print(f"Insertion count in BB_1:  {ins_count:,}")
        print(f"Deletion count in BB_1:   {del_count:,}")
        print(
            f"Reads trimmed below {args.min_adapter_start_pos}bp: {no_adapter_count:,}"
        )


if __name__ == "__main__":
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

    # Run params
    print(f"Adapter sequence:                   {args.adapter_seq}")
    print(f"Input fastq:                        {args.fq1_in}")
    print(f"Output fastq:                       {args.fq1_out}")
    print(f"Minimum adapter start position:     {args.min_adapter_start_pos}")
    print(f"Minimum adapter alignment score:    {args.min_align_score}")

    main(args)
