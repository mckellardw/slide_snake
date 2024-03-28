# Parallelized fastq splitting into N files, maintaining read order

# Requirements:
#   TODO

# Usage:
#   python splitNfqs.py in.fq.gz 12 12 False {optional: N_lines}

# TODO: add output options

# imports
import sys
import os

# import subprocess


# Function to extract out a section of the fastq file
def extractReads(fq_in, fq_out, start, stop):
    os.system(
        f"""
        zcat {fq_in} | sed -n '{start},{stop}p' > {fq_out}
        """
    )


# Count # of lines
def countLines(filename):
    with os.popen("zcat " + filename + "| wc -l") as f:
        return int(f.read().split()[0])


# Get arguments
fq = sys.argv[1]  # Path to fq file
n_chunks = int(sys.argv[2])  # Number of files to chunk the .fastq into
n_cores = int(sys.argv[3])  # Number of cores to parallelize with

try:  # Optional: number of total lines in fastq file (n_reads*4)
    n_total_lines = int(sys.argv[5])
except:
    n_total_lines = None

# Argument checks
## Make sure core count is lower/equal to chunk count
if n_cores > n_chunks:
    n_cores = n_chunks

## Count number of lines if not given
if n_total_lines is None:  # arg not given
    # n_total_lines = int(subprocess.check_output(['zcat' ,fq, '| wc -l']).split()[0])
    n_total_lines = countLines(fq)
elif n_total_lines // 4 != 0:  # Check that n_total_lines is divisible by 4
    n_total_lines = countLines(fq)

# Check line count...
if n_total_lines < 4:
    print(f"`{fq}` is empty! Exiting.")
    sys.exit()

# Compute start/stop pairs for reads
lines_per_file = n_total_lines // n_chunks
stops = [lines_per_file * n for n in range(1, n_chunks)]
stops = [
    stop - stop % 4 for stop in stops
]  # Ensure that stop positions aren't in the middle of a read (not divisible by 4)
stops.append(n_total_lines)  # Add final stop (EOF)

starts = [stops[i - 1] + 1 for i in range(1, n_chunks + 1)]  # Get start lines
starts.insert(0, 1)  # Add beginning of file to starts

if "fq" in fq:
    fq_out_list = [
        f"{fq.replace('.fq.gz','')}_{str(n).zfill(3)}.fq"
        for n in range(1, n_chunks + 1)
    ]
elif "fastq" in fq:
    fq_out_list = [
        f"{fq.replace('.fastq.gz','')}_{str(n).zfill(3)}.fastq"
        for n in range(1, n_chunks + 1)
    ]

# Extract reads and write new files
if n_cores > 1:  # Parallelize with `multiprocessing`
    import multiprocessing

    items = [(fq, fq_out_list[i], starts[i], stops[i]) for i in range(0, n_chunks)]
    with multiprocessing.Pool(n_cores) as pool:
        multi_out = pool.starmap(extractReads, items)
else:  # Single thread
    for i in list(range(0, n_chunks)):
        extractReads(fq, fq_out_list[i], starts[i], stops[i])
