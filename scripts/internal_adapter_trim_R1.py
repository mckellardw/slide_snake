# Script to remove the adapter sequence from R1

# Usage examples:
#   Snakemake:
#       python scripts/internal_adapter_trim_R1.py {params.INTERNAL_ADAPTER} {log} {threads} {params.TMPDIR} {input.MERGED_R1_FQ} {output.FINAL_R1_FQ}
#
#   bash:
#       python scripts/internal_adapter_trim_R1.py TCTTCAGCGTTCCCGAGA /workdir/dwm269/totalRNA/STRS-HD/data/align_out_rRNA/SH4/internal_adapter_trim_R1.log 18 /workdir/dwm269/totalRNA/STRS-HD/data/align_out_rRNA/SH4/tmp/seqtk /workdir/dwm269/totalRNA/STRS-HD/data/align_out_rRNA/SH4/tmp/SH4_R1_adapterTrim.fq.gz /workdir/dwm269/totalRNA/STRS-HD/data/align_out_rRNA/SH4/tmp/SH4_R1_finalInternalTrim.fq.gz

# imports
import sys
import os
import gzip
from Bio.SeqIO.QualityIO import FastqGeneralIterator
# from Bio.Blast import NCBIWWW
from Bio import pairwise2

# params
adapter_seq = sys.argv[1] # Curio's R1 adapter
log_out = sys.argv[2]  
n_cores = int(sys.argv[3])
tmp_dir = sys.argv[4] 
fq1_in = sys.argv[5]
fq1_out = sys.argv[6]

# fq2_in = sys.argv[4]
# fq2_out = sys.argv[5]

#TODO: add in read filtering for both R1 and R2 for reads which have alignment[0].score < 50

# Function to run internal trimming on a single .fq.gz file
def trim_fq(fq_in, fq_out, gz):
    # fq_in, fq_out, gz = args

    # Tallies for log file
    read_count = 0   # Read count
    ins_count = 0    # Insertion counter for BB_1
    del_count = 0    # Deletion counter for BB_1
    broken_count = 0 # Tally of reads that are removed (final read length <22bp)

    out_handle = open(fq_out, "w")

    if gz: 
        fq_iterator = FastqGeneralIterator(gzip.open(fq_in, "rt"))
    else:
        fq_iterator = FastqGeneralIterator(open(fq_in, "r"))

    for title, seq, qual in fq_iterator:
        read_count += 1

        ## match score = 4, mismatch = -0.5
        ## gap opening = -8, gap extension = -8
        alignment = pairwise2.align.localms(
            seq, 
            adapter_seq,
            4, -0.5, -8, -8,
            # penalize_end_gaps =[True, True],
            # score_only = True,
            one_alignment_only=True
        )

        #For troubleshooting:
        # if alignment[0].score > 55 and alignment[0].start > 9:
        #     print(alignment[0])
        #     print(pairwise2.format_alignment(*alignment[0]))
        #     print(seq[0:9]+seq[alignment[0].end:50])

        # Trim seq and qual
        start = alignment[0].start

        # Acount for reads with deletions in `BB_1`
        if start < 9:
            offset = 9-start
            seq_out = "N"*offset + seq[start:9] + seq[alignment[0].end:50]
            qual_out = "!"*offset + qual[start:9] + qual[alignment[0].end:50]

            del_count += 1
        else:
            if start > 9: 
                ins_count += 1

            ## Trim the first base in the read
            # seq_out = seq[start-8:start]+seq[alignment[0].end:50]
            # qual_out = qual[start-8:start]+qual[alignment[0].end:50]

            ## Trim the base closest to adapter
            seq_out = seq[0:9]+seq[alignment[0].end:50]
            qual_out = qual[0:9]+qual[alignment[0].end:50]

        if len(seq_out) < 22: broken_count += 1

        # Write to new .fq.gz file
        out_handle.write(
            f"@{title}\n{seq_out}\n+\n{qual_out}\n"
        )

    out_handle.close()

    return [read_count, ins_count, del_count, broken_count]

if n_cores > 1:
    # Source: https://superfastpython.com/multiprocessing-pool-for-loop/
    import multiprocessing
    # from joblib import Parallel, delayed

    # Split .fq file into {n_core} chunks
    os.system(
        f"""
        zcat {fq1_in} | seqkit split2 -p {n_cores} -O {tmp_dir} --force
        """
    )

    # Trim each chunked file, in parallel
    # tmp_fqs_in = [f"stdin.part_00{n}.fastq" for n in list(range(1,n_cores))]
    tmp_fqs_out = [f"{tmp_dir}/stdin.part_{str(n).zfill(3)}_trimmed.fastq" for n in list(range(1,n_cores))]

    items = [(f"{tmp_dir}/stdin.part_{str(n).zfill(3)}.fastq", f"{tmp_dir}/stdin.part_{str(n).zfill(3)}_trimmed.fastq", False) for n in list(range(1,n_cores))]
    print(items[1])
    with multiprocessing.Pool(n_cores) as pool:
        multi_out = pool.starmap(trim_fq, items)
    print("\n \n \n")
    print(multi_out)

    read_count = 0
    ins_count = 0
    del_count = 0
    broken_count = 0
    for i in list(range(0,len(multi_out))):
        read_count += multi_out[i][0]
        ins_count += multi_out[i][1]
        del_count += multi_out[i][2]
        broken_count += multi_out[i][3]

    # Merge and compress chunked/trimmed fastqs
    os.system(
        f"""
        cat {tmp_fqs_out} > {fq1_out.replace('.gz','')}
        pigz -p{n_cores} {fq1_out.replace('.gz','')}
        """
    )

    # Remove tmp fastqs
    for n in list(range(1,n_cores)):
        if os.path.isfile(f"stdin.part_{str(n).zfill(3)}.fastq"):
            os.remove()
        if os.path.isfile(f"stdin.part_{str(n).zfill(3)}_trimmed.fastq"):
            os.remove()

    with open(log_out, "w") as text_file:
        text_file.write(
            f"""
            Total read count:         {sum(read_count):,}
            Insertion count in BB_1:  {sum(ins_count):,}
            Deletion count in BB_2:   {sum(del_count):,}
            Reads trimmed below 22bp: {sum(broken_count):,}
            """
        )
else:
    out = trim_fq(fq1_in, fq1_out, True)

    read_count, ins_count, del_count, broken_count = out

    # Compress trimmed fastq
    os.system(
        f"""
        pigz -p{n_cores} {fq1_out.replace('.gz','')}
        """
    )

    with open(log_out, "w") as text_file:
        text_file.write(
            f"""
            Total read count:         {read_count:,}
            Insertion count in BB_1:  {ins_count:,}
            Deletion count in BB_2:   {del_count:,}
            Reads trimmed below 22bp: {broken_count:,}
            """
        )

print("Done!")