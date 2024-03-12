#!/bin/bash

# Function to display usage information
usage() {
    echo "Usage: $0 -r1 <R1_FQ_IN> -r2 <R2_FQ_IN> -o1 <R1_FQ_OUT> -o2 <R2_FQ_OUT> -p <BC_PATTERN> -w <BC_WHITELIST> -l <LOG>"
    echo
    echo "Options:"
    echo " --r1 <R1_FQ_IN>     Input file for R1 reads"
    echo " --r2 <R2_FQ_IN>     Input file for R2 reads"
    echo " --o1 <R1_FQ_OUT>    Output file for R1 reads"
    echo " --o2 <R2_FQ_OUT>    Output file for R2 reads"
    echo " -p <BC_PATTERN>    Barcode pattern"
    echo " -w <BC_WHITELIST>  Barcode whitelist file"
    echo " -l <LOG>           Log file"
    exit 1
}
# Initialize variables
R1_FQ_IN=""
R2_FQ_IN=""
R1_FQ_OUT=""
R2_FQ_OUT=""
BC_PATTERN=""
BC_WHITELIST=""
LOG=""

# Parse arguments
while (( "$#" )); do
    case "$1" in
        -r1)
            R1_FQ_IN="$2"
            shift 2
            ;;
        -r2)
            R2_FQ_IN="$2"
            shift 2
            ;;
        -o1)
            R1_FQ_OUT="$2"
            shift 2
            ;;
        -o2)
            R2_FQ_OUT="$2"
            shift 2
            ;;
        -p)
            BC_PATTERN="$2"
            shift 2
            ;;
        -w)
            BC_WHITELIST="$2"
            shift 2
            ;;
        -l)
            LOG="$2"
            shift 2
            ;;
        *)
            echo "Error: Unknown option $1"
            usage
            ;;
    esac
done

# Check if all required arguments are provided
# if [ -z "$R1_FQ_IN" ] || [ -z "$R2_FQ_IN" ] || [ -z "$R1_FQ_OUT" ] || [ -z "$R2_FQ_OUT" ] || [ -z "$BC_PATTERN" ] || [ -z "$BC_WHITELIST" ] || [ -z "$LOG" ]; then
#     echo "Error: Missing required arguments."
#     usage
# fi

# Check if all required arguments are provided
if [ -z "$R1_FQ_IN" ]; then
    echo "Error: Missing argument for R1 input file."
    usage
fi
if [ -z "$R2_FQ_IN" ]; then
    echo "Error: Missing argument for R2 input file."
    usage
fi
if [ -z "$R1_FQ_OUT" ]; then
    echo "Error: Missing argument for R1 output file."
    usage
fi
if [ -z "$R2_FQ_OUT" ]; then
    echo "Error: Missing argument for R2 output file."
    usage
fi
if [ -z "$BC_PATTERN" ]; then
    echo "Error: Missing argument for barcode pattern."
    usage
fi
if [ -z "$BC_WHITELIST" ]; then
    echo "Error: Missing argument for barcode whitelist file."
    usage
fi
if [ -z "$LOG" ]; then
    echo "Error: Missing argument for log file."
    usage
fi

# umitools extract | https://umi-tools.readthedocs.io/en/latest/reference/extract.html
umi_tools extract \
    --extract-method=string \
    --bc-pattern=${BC_PATTERN} \
    --whitelist=${BC_WHITELIST} \
    --stdin=${R1_FQ_IN} \
    --read2-in=${R2_FQ_IN} \
    --stdout=${R1_FQ_OUT} \
    --read2-out=${R2_FQ_OUT} \
    --log=${LOG}
