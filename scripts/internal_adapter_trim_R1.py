# Script to remove the adapter sequence from R1

# Usage examples:
#   Snakemake:
#       python scripts/internal_adapter_trim_R1.py {params.INTERNAL_ADAPTER} {log} {threads} {params.TMPDIR} {input.MERGED_R1_FQ} {output.FINAL_R1_FQ}
#
#   bash:
#       python scripts/internal_adapter_trim_R1.py TCTTCAGCGTTCCCGAGA /workdir/dwm269/totalRNA/STRS-HD/data/align_out_rRNA/SH4/internal_adapter_trim_R1.log 18 /workdir/dwm269/totalRNA/STRS-HD/data/align_out_rRNA/SH4/tmp/seqkit /workdir/dwm269/totalRNA/STRS-HD/data/align_out_rRNA/SH4/tmp/SH4_R1_adapterTrim.fq.gz /workdir/dwm269/totalRNA/STRS-HD/data/align_out_rRNA/SH4/tmp/SH4_R1_finalInternalTrim.fq.gz

# R1 structure:
#  [BB1- 8bp][Adapter Seq - 18bp][BB2|6bp][UMI-7bp][TTTTTTTTTTT]

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
# prefix = sys.argv[5]

# fq2_in = sys.argv[7]
# fq2_out = sys.argv[8]

#TODO: add in read filtering for both R1 and R2 for reads which have alignment[0].score < 50
#TODO: fix indentation in log files

#TODO: add fq1_in param check for fq vs fastq suffix; replace all '.fq. with `suffix`

# Function to run internal trimming on a single .fq.gz file
#TODO- update to newer biopython functions...
# warnings.warn(
# /home/dwm269/miniconda3/envs/STARsolo/lib/python3.9/site-packages/Bio/pairwise2.py:278: BiopythonDeprecationWarning: Bio.pairwise2 
# has been deprecated, and we intend to remove it in a future release of Biopython. As an alternative, please consider using 
# Bio.Align.PairwiseAligner as a replacement, and contact the Biopython developers if you still need the Bio.pairwise2 module.
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

        # Useful resource: https://www.bioinformaticscrashcourse.com/10.1_Alignment.html
        ## match score = 4, mismatch = -0.5
        ## gap opening = -6, gap extension = -6
        alignment = pairwise2.align.localms(
            seq, 
            adapter_seq,
            4, -0.5, -6, -6,
            # score_only = True,
            one_alignment_only=True
        )

        #For troubleshooting:
        # if alignment[0].score < 60 and alignment[0].score > 55:# and alignment[0].start > 9:
        #     print(alignment[0])
        #     print(pairwise2.format_alignment(*alignment[0]))
        #     print(seq[0:8]+seq[alignment[0].end:])

        # Trim seq and qual
        start = alignment[0].start

        # Acount for reads with deletions in `BB_1`
        if start < 8: # Deletion in BB_1
            offset = 8-start
            seq_out = "N"*offset + seq[start:8] + seq[alignment[0].end:]
            qual_out = "!"*offset + qual[start:8] + qual[alignment[0].end:]

            del_count += 1
        else:
            if start > 8: # Insertion in BB_1
                ins_count += 1

            ## Trim the first base in the read
            # seq_out = seq[start-8:start]+seq[alignment[0].end:]
            # qual_out = qual[start-8:start]+qual[alignment[0].end:]

            ## Trim the base closest to adapter
            seq_out = seq[0:8]+seq[alignment[0].end:]
            qual_out = qual[0:8]+qual[alignment[0].end:]

        # Broken read; erase R1 and add `NNNNNNNNNN` with qval=0 (!!!!!!!!!!)
        # Alignment score cuttoff was determined by spot-checking ~100 bead barcodes, based on predicted adapter location
        if len(seq_out) < 22 or alignment[0].score < 58: 
            seq_out = "N"*10 
            qual_out = "!"*10
            broken_count += 1

        # Write to new .fq.gz file
        out_handle.write(
            f"@{title}\n{seq_out}\n+\n{qual_out}\n"
        )

    out_handle.close()

    return [read_count, ins_count, del_count, broken_count]

if n_cores > 1:
    # Source: https://superfastpython.com/multiprocessing-pool-for-loop/
    import multiprocessing

    # Split .fq file into {n_core} chunks
    os.system(
        # Super high mem usage
        # f""" 
        # zcat {fq1_in} | seqkit split -p {n_cores} -O {tmp_dir} --force
        # """
        
        # Custom script to split .fq w/ sed
        f""" 
        python scripts/splitNfqs.py {fq1_in} {n_cores} {n_cores} "False"
        """
    )

    # Trim each chunked file, in parallel
    # tmp_fqs_in = [f"stdin.part_00{n}.fastq" for n in list(range(1,n_cores+1))] #seqkit split
    # tmp_fqs_in = [f"{fq1_in.replace('fq.gz','')}_{str(n).zfill(3)}.fq" for n in range(1,n_chunks+1)] #splitNfqs

    tmp_fqs_out = [f"{fq1_in.replace('.fq.gz','')}_{str(n).zfill(3)}_trimmed.fq" for n in list(range(1,n_cores+1))]
    items = [(f"{fq1_in.replace('.fq.gz','')}_{str(n).zfill(3)}.fq", f"{fq1_in.replace('.fq.gz','')}_{str(n).zfill(3)}_trimmed.fq", False) for n in list(range(1,n_cores+1))]
    
    with multiprocessing.Pool(n_cores) as pool:
        multi_out = pool.starmap(trim_fq, items)

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
        cat {fq1_in.replace('.fq.gz','')}_*_trimmed.fq > {fq1_out.replace('.gz','')}
        pigz -p{n_cores} {fq1_out.replace('.gz','')}
        """
    )
    # {tmp_dir}/stdin.part_*_trimmed.fastq

    # Remove tmp fastqs
    # for n in list(range(1,n_cores+1)):
    #     if os.path.isfile(f"{tmp_dir}/stdin.part_{str(n).zfill(3)}.fastq"):
    #         os.remove()
    #     if os.path.isfile(f"{tmp_dir}/stdin.part_{str(n).zfill(3)}_trimmed.fastq"):
    #         os.remove()

    if os.path.isfile(fq1_out):
        os.system(
            f"""
            rm {fq1_in.replace('.fq.gz','')}_[0-9][0-9][0-9].fq {fq1_in.replace('.fq.gz','')}_*_trimmed.fq
            """
        )
    # {tmp_dir}/stdin.part_*_trimmed.fastq ||:

    with open(log_out, "w") as text_file:
        text_file.write(
            f"""
Total read count:         {read_count:,}
Insertion count in BB_1:  {ins_count:,}
Deletion count in BB_2:   {del_count:,}
Reads trimmed below 22bp: {broken_count:,}
            """
        )
else:
    out = trim_fq(fq1_in, fq1_out, True)

    read_count, ins_count, del_count, broken_count = out

    # Compress trimmed fastq
    os.system(
        f"""
        pigz -f -p{n_cores} {fq1_out.replace('.gz','')}
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