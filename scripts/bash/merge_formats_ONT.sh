#!/bin/bash

set -eo pipefail

# Function to display help message
usage() {
    cat << EOF
Usage: $0 -d TMPDIR -r ONT_reads -c CHUNK_SIZE -o MERGED_FQ -l LOG -t THREADS [-f output_format]
       -d TMPDIR        : Temporary directory
       -r ONT_reads     : Space-delimited list of ONT reads or a regex pattern
       -c CHUNK_SIZE    : Chunk size
       -o MERGED_FQ     : Output merged file (with .fq.gz/.fastq.gz for FASTQ or .bam for BAM)
       -l LOG           : Log file
       -t THREADS       : Number of threads
       -f output_format : Output format (fastq or bam, default: fastq)
EOF
    exit 1
}

# Function to log messages
log_message() {
    echo "$(date '+%Y-%m-%d %H:%M:%S') - $1" | tee -a "$LOG"
}

# Function to check if an argument is missing
check_arg() {
    if [ -z "$2" ]; then
        echo "Error: Missing required argument -$1 $3"
        return 1
    fi
    return 0
}


# Default output format
OUTPUT_FORMAT="fastq"

while getopts ":d:r:c:o:l:t:f:" opt; do
    case $opt in
        d) TMPDIR=$OPTARG ;;
        r) ONT_reads=$OPTARG ;;
        c) CHUNK_SIZE=$OPTARG ;;
        o) MERGED_FQ=$OPTARG ;;
        l) LOG=$OPTARG ;;
        t) THREADS=$OPTARG ;;
        f) OUTPUT_FORMAT=$OPTARG ;;
        \?) echo "Invalid option: -$OPTARG" >&2; usage ;;
        :) echo "Option -$OPTARG requires an argument." >&2; usage ;;
    esac
done

# Check for missing arguments
missing_args=0
check_arg "d" "$TMPDIR" "TMPDIR" || missing_args=$((missing_args + 1))
check_arg "r" "$ONT_reads" "ONT_reads" || missing_args=$((missing_args + 1))
check_arg "c" "$CHUNK_SIZE" "CHUNK_SIZE" || missing_args=$((missing_args + 1))
check_arg "o" "$MERGED_FQ" "MERGED_FQ" || missing_args=$((missing_args + 1))
check_arg "l" "$LOG" "LOG" || missing_args=$((missing_args + 1))
check_arg "t" "$THREADS" "THREADS" || missing_args=$((missing_args + 1))

if [ $missing_args -gt 0 ]; then
    echo "Total missing arguments: $missing_args"
    usage
    exit 1
fi


# Validate output format
if [ "$OUTPUT_FORMAT" != "fastq" ] && [ "$OUTPUT_FORMAT" != "bam" ]; then
    echo "Error: Invalid output format. Must be 'fastq' or 'bam'."
    exit 1
fi

# Create TMPDIR and LOG file
if ! mkdir -p "$TMPDIR"; then
    echo "Error: Unable to create temporary directory $TMPDIR"
    exit 1
fi

if ! touch "$LOG"; then
    echo "Error: Unable to create log file $LOG"
    exit 1
fi

mkdir -p $TMPDIR
# rm -rf $TMPDIR/* # clear out tmp dir if needed

mkdir -p $(dirname $LOG)
: > $LOG

printf 'TMPDIR:           %s\n' "${TMPDIR}">> ${LOG}
echo "ONT_reads:        ${ONT_reads}" >> ${LOG}
echo "CHUNK_SIZE:       ${CHUNK_SIZE}" >> ${LOG}
echo "MERGED_FQ:        ${MERGED_FQ}" >> ${LOG}
echo "LOG:              ${LOG}" >> ${LOG}
echo "THREADS:          ${THREADS}" >> ${LOG}
echo "OUTPUT_FORMAT:    ${OUTPUT_FORMAT}" >> ${LOG}
echo "" >> $LOG  # Add a blank line for readability


log_message "Read files:" 
for f in $ONT_reads; do
    log_message "   $f" 
done

F_base=${MERGED_FQ%.gz}
F_base=${F_base%.bam}

if [[ $(echo $ONT_reads | wc -w) -eq 1 && ! $ONT_reads =~ "*" ]]; then
    F=$ONT_reads
    if [[ $OUTPUT_FORMAT == "fastq" ]]; then
        if [[ $F == *.fq.gz || $F == *.fastq.gz ]]; then
            cp $F $MERGED_FQ 2>> $LOG
        elif [[ $F == *.sam || $F == *.bam || $F == *.cram ]]; then
            if [[ -f $F ]]; then
                samtools fastq $F > $F_base 2>> $LOG
            else
                echo "File [ $F ] does not exist." >> $LOG
            fi
            pigz -p$THREADS $F_base 2>> $LOG
        else
            echo "File type for [$F] not supported!" >> $LOG
        fi
    elif [[ $OUTPUT_FORMAT == "bam" ]]; then
        if [[ $F == *.bam || $F == *.cram ]]; then
            cp $F $MERGED_FQ 2>> $LOG
        elif [[ $F == *.sam ]]; then
            samtools view -bS $F > $MERGED_FQ 2>> $LOG
        else
            echo "File type for [$F] not supported!" >> $LOG
        fi
    else
        echo "Output format [$OUTPUT_FORMAT] not supported!" >> $LOG
    fi
elif [[ $(echo $ONT_reads | wc -w) -eq 1 && $ONT_reads =~ "*" ]]; then
    F_list=($(ls $ONT_reads))
    echo 'Regex-ed file list:' >> $LOG
    for f in "${F_list[@]}"; do
        echo "   $f" >> $LOG
    done

    if [[ $OUTPUT_FORMAT == "fastq" ]]; then
        for i in "${!F_list[@]}"; do
            F=${F_list[$i]}
            if [[ $F == *.fq.gz || $F == *.fastq.gz ]]; then
                echo "Adding $F to output fastq" >> $LOG
                zcat $F >> $F_base
            elif [[ $F == *.sam || $F == *.bam || $F == *.cram ]]; then
                samtools fastq $F >> $F_base 2>> $LOG
            else
                echo "File type for [$F] not supported!" >> $LOG
            fi
        done
        pigz -p$THREADS $F_base
    elif [[ $OUTPUT_FORMAT == "bam" ]]; then
        for i in "${!F_list[@]}"; do
            F=${F_list[$i]}
            if [[ $F == *.bam || $F == *.cram ]]; then
                echo "Adding $F to output bam" >> $LOG
                samtools cat -o $F_base.bam $F_base.bam $F 2>> $LOG
            elif [[ $F == *.sam ]]; then
                samtools view -bS $F >> $F_base.bam 2>> $LOG
            else
                echo "File type for [$F] not supported!" >> $LOG
            fi
        done
        mv $F_base.bam $MERGED_FQ
    else
        echo "Output format [$OUTPUT_FORMAT] not supported!" >> $LOG
    fi
else
    if [[ $OUTPUT_FORMAT == "fastq" ]]; then
        for F in $ONT_reads; do
            if [[ $F == *.fq.gz || $F == *.fastq.gz ]]; then
                if [[ -f $F ]]; then
                    echo "Adding $F to output fastq" >> $LOG
                    zcat $F >> $F_base
                else
                    echo "File [ $F ] does not exist." >> $LOG
                fi
            elif [[ $F == *.sam || $F == *.bam || $F == *.cram ]]; then
                if [[ -f $F ]]; then
                    samtools fastq $F >> $F_base 2>> $LOG
                else
                    echo "File [ $F ] does not exist." >> $LOG
                fi
            else
                echo "File type for [$F] not supported!" >> $LOG
            fi
        done
        pigz -p$THREADS $F_base 2>> $LOG
    elif [[ $OUTPUT_FORMAT == "bam" ]]; then
        for F in $ONT_reads; do
            if [[ $F == *.bam || $F == *.cram ]]; then
                if [[ -f $F ]]; then
                    echo "Adding $F to output bam" >> $LOG
                    samtools cat -o $F_base.bam $F_base.bam $F 2>> $LOG
                elif [[ $F == *.sam ]]; then
                    samtools view -bS $F >> $F_base.bam 2>> $LOG
                else
                    echo "File type for [$F] not supported!" >> $LOG
                fi
            else
                echo "File type for [$F] not supported!" >> $LOG
            fi
        done
        mv $F_base.bam $MERGED_FQ
    else
        echo "Output format [$OUTPUT_FORMAT] not supported!" >> $LOG
    fi
fi

log_message "Done!"