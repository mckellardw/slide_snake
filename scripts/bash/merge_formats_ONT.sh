#!/bin/bash

set -euo pipefail

# Function to display help message
usage() {
    cat << EOF
Usage: $0 -d TMPDIR -r ONT_reads -c CHUNK_SIZE -o MERGED_FQ -t THREADS [-f output_format]
       -d TMPDIR        : Temporary directory
       -r ONT_reads     : Space-delimited list of ONT reads or a regex pattern
       -c CHUNK_SIZE    : Chunk size
       -o MERGED_FQ     : Output merged file (with .fq.gz/.fastq.gz for FASTQ or .bam for BAM)
       -t THREADS       : Number of threads
       -f output_format : Output format (fastq or bam, default: fastq)
EOF
    exit 1
}

# Function to log messages
log_message() {
    local message="$1"
    printf "%s - %s\n" "$(date '+%Y-%m-%d %H:%M:%S')" "$message"
}

# Function to check if an argument is missing
check_arg() {
    local arg_name="$1"
    local arg_value="$2"
    local arg_desc="$3"
    if [ -z "$arg_value" ]; then
        printf "Error: Missing required argument -%s %s\n" "$arg_name" "$arg_desc"
        return 1
    fi
    return 0
}

# Default output format
OUTPUT_FORMAT="fastq"

while getopts ":d:r:c:o:t:f:" opt; do
    case $opt in
        d) TMPDIR=$OPTARG ;;
        r) ONT_reads=$OPTARG ;;
        c) CHUNK_SIZE=$OPTARG ;;
        o) MERGED_FQ=$OPTARG ;;
        t) THREADS=$OPTARG ;;
        f) OUTPUT_FORMAT=$OPTARG ;;
        \?) printf "Invalid option: -%s\n" "$OPTARG" >&2; usage ;;
        :) printf "Option -%s requires an argument.\n" "$OPTARG" >&2; usage ;;
    esac
done

# Check for missing arguments
missing_args=0
check_arg "d" "$TMPDIR" "TMPDIR" || missing_args=$((missing_args + 1))
check_arg "r" "$ONT_reads" "ONT_reads" || missing_args=$((missing_args + 1))
check_arg "c" "$CHUNK_SIZE" "CHUNK_SIZE" || missing_args=$((missing_args + 1))
check_arg "o" "$MERGED_FQ" "MERGED_FQ" || missing_args=$((missing_args + 1))
check_arg "t" "$THREADS" "THREADS" || missing_args=$((missing_args + 1))

if [ $missing_args -gt 0 ]; then
    printf "Total missing arguments: %d\n" "$missing_args"
    usage
    exit 1
fi

# Validate output format
if [ "$OUTPUT_FORMAT" != "fastq" ] && [ "$OUTPUT_FORMAT" != "bam" ]; then
    printf "Error: Invalid output format. Must be 'fastq' or 'bam'.\n"
    exit 1
fi

# Create TMPDIR
if ! mkdir -p "$TMPDIR"; then
    printf "Error: Unable to create temporary directory %s\n" "$TMPDIR"
    exit 1
fi

printf 'TMPDIR:           %s\n' "${TMPDIR}"
printf 'ONT_reads:        %s\n' "${ONT_reads}"
printf 'CHUNK_SIZE:       %s\n' "${CHUNK_SIZE}"
printf 'MERGED_FQ:        %s\n' "${MERGED_FQ}"
printf 'THREADS:          %s\n' "${THREADS}"
printf 'OUTPUT_FORMAT:    %s\n' "${OUTPUT_FORMAT}"
printf '\n'  # Add a blank line for readability

log_message "Read files:" 
for f in $ONT_reads; do
    log_message "   $f" 
done

F_base=${MERGED_FQ%.gz}
F_base=${F_base%.bam}

if [[ $(echo $ONT_reads | wc -w) -eq 1 && ! $ONT_reads =~ "*" ]]; then
    F=$ONT_reads
    case $OUTPUT_FORMAT in
        fastq)
            case $F in
                *.fq.gz|*.fastq.gz)
                    cp "$F" "$MERGED_FQ"
                    ;;
                *.sam|*.bam|*.cram)
                    if [[ -f $F ]]; then
                        samtools fastq "$F" > "$F_base"
                    else
                        printf "File [ %s ] does not exist.\n" "$F"
                    fi
                    pigz -p"$THREADS" "$F_base"
                    ;;
                *)
                    printf "File type for [%s] not supported!\n" "$F"
                    ;;
            esac
            ;;
        bam)
            case $F in
                *.bam|*.cram)
                    cp "$F" "$MERGED_FQ"
                    ;;
                *.sam)
                    samtools view -bS "$F" > "$MERGED_FQ"
                    ;;
                *)
                    printf "File type for [%s] not supported!\n" "$F"
                    ;;
            esac
            ;;
        *)
            printf "Output format [%s] not supported!\n" "$OUTPUT_FORMAT"
            ;;
    esac
elif [[ $(echo $ONT_reads | wc -w) -eq 1 && $ONT_reads =~ "*" ]]; then
    F_list=($(ls $ONT_reads))
    printf 'Regex-ed file list:\n'
    for f in "${F_list[@]}"; do
        printf "   %s\n" "$f"
    done

    case $OUTPUT_FORMAT in
        fastq)
            for F in "${F_list[@]}"; do
                case $F in
                    *.fq.gz|*.fastq.gz)
                        printf "Adding %s to output fastq\n" "$F"
                        zcat "$F" >> "$F_base"
                        ;;
                    *.sam|*.bam|*.cram)
                        samtools fastq "$F" >> "$F_base"
                        ;;
                    *)
                        printf "File type for [%s] not supported!\n" "$F"
                        ;;
                esac
            done
            pigz -p"$THREADS" "$F_base"
            ;;
        bam)
            for F in "${F_list[@]}"; do
                case $F in
                    *.bam|*.cram)
                        printf "Adding %s to output bam\n" "$F"
                        samtools cat -o "$F_base.bam" "$F_base.bam" "$F"
                        ;;
                    *.sam)
                        samtools view -bS "$F" >> "$F_base.bam"
                        ;;
                    *)
                        printf "File type for [%s] not supported!\n" "$F"
                        ;;
                esac
            done
            mv "$F_base.bam" "$MERGED_FQ"
            ;;
        *)
            printf "Output format [%s] not supported!\n" "$OUTPUT_FORMAT"
            ;;
    esac
else
    case $OUTPUT_FORMAT in
        fastq)
            for F in $ONT_reads; do
                case $F in
                    *.fq.gz|*.fastq.gz)
                        if [[ -f $F ]]; then
                            printf "Adding %s to output fastq\n" "$F"
                            zcat "$F" >> "$F_base"
                        else
                            printf "File [ %s ] does not exist.\n" "$F"
                        fi
                        ;;
                    *.sam|*.bam|*.cram)
                        if [[ -f $F ]]; then
                            samtools fastq "$F" >> "$F_base"
                        else
                            printf "File [ %s ] does not exist.\n" "$F"
                        fi
                        ;;
                    *)
                        printf "File type for [%s] not supported!\n" "$F"
                        ;;
                esac
            done
            pigz -p"$THREADS" "$F_base"
            ;;
        bam)
            for F in $ONT_reads; do
                case $F in
                    *.bam|*.cram)
                        if [[ -f $F ]]; then
                            printf "Adding %s to output bam\n" "$F"
                            samtools cat -o "$F_base.bam" "$F_base.bam" "$F"
                        elif [[ $F == *.sam ]]; then
                            samtools view -bS "$F" >> "$F_base.bam"
                        else
                            printf "File type for [%s] not supported!\n" "$F"
                        fi
                        ;;
                    *)
                        printf "File type for [%s] not supported!\n" "$F"
                        ;;
                esac
            done
            mv "$F_base.bam" "$MERGED_FQ"
            ;;
        *)
            printf "Output format [%s] not supported!\n" "$OUTPUT_FORMAT"
            ;;
    esac
fi

log_message "Done!"