#!/bin/bash

set -euo pipefail

# Function to display help message
usage() {
    cat << EOF
Usage: $0 -d TMPDIR -r READS -o MERGED_R1_FQ -p MERGED_R2_FQ -t THREADS
       -d TMPDIR        : Temporary directory
       -r READS         : Space-delimited list of paired-end reads or a regex pattern
       -o MERGED_R1_FQ  : Output merged R1 file
       -p MERGED_R2_FQ  : Output merged R2 file
       -t THREADS       : Number of threads
EOF
    exit 1
}

# Function to log messages
log_message() {
    local message="$1"
    printf "%s - %s\n" "$(date '+%Y-%m-%d %H:%M:%S')" "$message"
}

# Parse arguments
while getopts ":d:r:o:p:t:" opt; do
    case $opt in
        d) TMPDIR=$OPTARG ;;
        r) READS=$OPTARG ;;
        o) MERGED_R1_FQ=$OPTARG ;;
        p) MERGED_R2_FQ=$OPTARG ;;
        t) THREADS=$OPTARG ;;
        \?) printf "Invalid option: -%s\n" "$OPTARG" >&2; usage ;;
        :) printf "Option -%s requires an argument.\n" "$OPTARG" >&2; usage ;;
    esac
done

# Check for missing arguments
if [ -z "${TMPDIR:-}" ] || [ -z "${READS:-}" ] || [ -z "${MERGED_R1_FQ:-}" ] || [ -z "${MERGED_R2_FQ:-}" ] || [ -z "${THREADS:-}" ]; then
    usage
fi

# Create TMPDIR
if ! mkdir -p "$TMPDIR"; then
    printf "Error: Unable to create temporary directory %s\n" "$TMPDIR"
    exit 1
fi

log_message "Temporary directory: $TMPDIR"
log_message "Reads: $READS"
log_message "Output R1: $MERGED_R1_FQ"
log_message "Output R2: $MERGED_R2_FQ"
log_message "Threads: $THREADS"

# Handle regex or space-delimited input
if [[ $READS =~ "*" ]]; then
    R1_LIST=($(ls $READS | grep "_R1"))
    R2_LIST=($(ls $READS | grep "_R2"))
else
    R1_LIST=($(echo $READS | tr ' ' '\n' | grep "_R1"))
    R2_LIST=($(echo $READS | tr ' ' '\n' | grep "_R2"))
fi

# Validate paired-end files
if [ ${#R1_LIST[@]} -ne ${#R2_LIST[@]} ]; then
    printf "Error: Mismatched number of R1 and R2 files.\n"
    exit 1
fi

log_message "Merging R1 files..."
zcat "${R1_LIST[@]}" > "${MERGED_R1_FQ%.gz}"
pigz -p "$THREADS" "${MERGED_R1_FQ%.gz}"

log_message "Merging R2 files..."
zcat "${R2_LIST[@]}" > "${MERGED_R2_FQ%.gz}"
pigz -p "$THREADS" "${MERGED_R2_FQ%.gz}"

log_message "Merging completed successfully!"
