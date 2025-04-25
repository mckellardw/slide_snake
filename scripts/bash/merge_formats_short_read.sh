#!/bin/bash

set -euo pipefail

# Function to display help message
usage() {
    cat << EOF
Usage: $0 -d TMPDIR -1 R1_LIST -2 R2_LIST -o MERGED_R1_FQ -p MERGED_R2_FQ -t THREADS
       -d TMPDIR        : Temporary directory
       -1 R1_LIST       : Space-delimited list of R1 files
       -2 R2_LIST       : Space-delimited list of R2 files
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
while getopts ":d:1:2:o:p:t:" opt; do
    case $opt in
        d) TMPDIR=$OPTARG ;;
        1) R1_LIST=$OPTARG ;;
        2) R2_LIST=$OPTARG ;;
        o) MERGED_R1_FQ=$OPTARG ;;
        p) MERGED_R2_FQ=$OPTARG ;;
        t) THREADS=$OPTARG ;;
        \?) printf "Invalid option: -%s\n" "$OPTARG" >&2; usage ;;
        :) printf "Option -%s requires an argument.\n" "$OPTARG" >&2; usage ;;
    esac
done

# Check for missing arguments
if [ -z "${TMPDIR:-}" ] || [ -z "${R1_LIST:-}" ] || [ -z "${R2_LIST:-}" ] || [ -z "${MERGED_R1_FQ:-}" ] || [ -z "${MERGED_R2_FQ:-}" ] || [ -z "${THREADS:-}" ]; then
    usage
fi

# Create TMPDIR
if ! mkdir -p "$TMPDIR"; then
    printf "Error: Unable to create temporary directory %s\n" "$TMPDIR"
    exit 1
fi

log_message "Temporary directory: $TMPDIR"
log_message "R1 files: $R1_LIST"
log_message "R2 files: $R2_LIST"
log_message "Output R1: $MERGED_R1_FQ"
log_message "Output R2: $MERGED_R2_FQ"
log_message "Threads: $THREADS"

# Convert space-delimited lists into arrays
R1_ARRAY=($R1_LIST)
R2_ARRAY=($R2_LIST)

# Validate paired-end files
if [ ${#R1_ARRAY[@]} -ne ${#R2_ARRAY[@]} ]; then
    printf "Error: Mismatched number of R1 and R2 files.\n"
    exit 1
fi

# Handle R1 files
if [ ${#R1_ARRAY[@]} -eq 1 ]; then
    log_message "Only one R1 file provided. Copying directly to output."
    cp "${R1_ARRAY[0]}" "$MERGED_R1_FQ"
else
    log_message "Merging R1 files..."
    zcat "${R1_ARRAY[@]}" > "${MERGED_R1_FQ%.gz}"
    pigz -p "$THREADS" "${MERGED_R1_FQ%.gz}"
fi

# Handle R2 files
if [ ${#R2_ARRAY[@]} -eq 1 ]; then
    log_message "Only one R2 file provided. Copying directly to output."
    cp "${R2_ARRAY[0]}" "$MERGED_R2_FQ"
else
    log_message "Merging R2 files..."
    zcat "${R2_ARRAY[@]}" > "${MERGED_R2_FQ%.gz}"
    pigz -p "$THREADS" "${MERGED_R2_FQ%.gz}"
fi

log_message "Merging completed successfully!"
