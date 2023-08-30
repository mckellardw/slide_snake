#!/bin/bash

# Example usage:
# bash script.sh mature.fa "Mus musculus" mmu mature 2> mmu/bt2_mature.log


# Input parameters
FASTA_FILE=$1
SPECIES=$2
OUTDIR=$3
PREFIX=$4

mkdir -p $OUTDIR

# Replace spaces in species name with underscores
SPECIES_NO_SPACES=$(echo $SPECIES | tr ' ' '_')

# Filter fasta file for specific species, and change Us to Ts for DNA seq data
cat $FASTA_FILE \
| awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' \
| grep -A 1 "$SPECIES" \
| sed '/^[^>]/s/U/T/g' \
| awk '/^>/{sub(/.*/, ">"$5)}1' \
> $OUTDIR/${SPECIES_NO_SPACES}_${PREFIX}_DNA.fa

# Generate bowtie2 index
bowtie2-build $OUTDIR/${SPECIES_NO_SPACES}_${PREFIX}_DNA.fa $OUTDIR/$PREFIX
