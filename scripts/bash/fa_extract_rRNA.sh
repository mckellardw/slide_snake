#!/bin/bash

# Check if the correct number of arguments was provided
if [ "$#" -ne 2 ]; then
    echo "Usage: scripts/bash/extract_rRNA_fasta.sh FASTA_cDNA FASTA_rRNA"
    exit 1
fi

# Assign the first argument to FASTA_cDNA and the second to FASTA_rRNA
FASTA_cDNA=$1
FASTA_rRNA=$2

# Fields for fasta are different for GENCODE
zcat ${FASTA_cDNA} \
| awk \
    -v RS="\n>" \
    -v FS=" " '$5=="gene_biotype:rRNA" || $5=="gene_biotype:Mt_rRNA" { print ">"$0 }' \
| gzip \
> ${FASTA_rRNA}