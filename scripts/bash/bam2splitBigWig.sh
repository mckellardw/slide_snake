#!/bin/bash

# Usage:
# bash bam2splitBigWig /path/to/file.bam num_cores /path/to/outdir /path/to/genome.chrom.sizes

# bam2splitBigWig (v1.0) - input a .bam file, output two bigWig files (one for each strand, positive/negative)
#                        - Also outputs bigwigs in log scale (IGV can't handle log w/ neg strand)

INBAM=$1 # path to .bam file (sorted & indexed already!)
CORE=$2 # number of cores for parallelization
OUTDIR=$3 # output directory
CHINFO=$4 # path to genome chrom.sizes

#TODO: add flags/options for including merged outputs and log10 outputs

PREFIX=`echo ${INBAM} | rev | cut -d / -f 1 | cut -d . -f 2- | rev`

echo ".bam file location:    "${INBAM}
echo "Max cores:             "${CORE}
echo "Output location:       "${OUTDIR}
echo "Chromosome size info:  "${CHINFO}
echo
echo

mkdir -p ${OUTDIR}
touch ${OUTDIR}/bam2splitBigWig.kill.warnings

# convert bam to bed &  split for strandedness
echo "Converting to bed and splitting strands..."
# bedtools bamtobed -i ${INBAM} | gzip > ${OUTDIR}/${PREFIX}.bed.gz

bedtools bamtobed -i ${INBAM}  2> ${OUTDIR}/bam2splitBigWig.kill.warnings | \
  awk 'BEGIN{OFS="\t"} ($5 > 0){print $0}' | \
  awk 'BEGIN{OFS="\t"} ($6 == "+") {print $1,$2,$2+1,$4,$5,$6}; ($6 == "-") {print $1,$3-1,$3,$4,$5,$6}' | \
  gzip > ${OUTDIR}/${PREFIX}.bed.gz

## Convert to bedGraph ... Cannot gzip these, unfortunately.

# echo "Converting to bedGraph...\n"
# Not stranded
# bedtools genomecov -bg -i ${OUTDIR}/${PREFIX}.bed.gz -g ${CHINFO} > ${OUTDIR}/${PREFIX}.bedGraph

echo "Converting to bedGraph..."
# Stranded
bedtools genomecov -bg -i ${OUTDIR}/${PREFIX}.bed.gz -g ${CHINFO} -strand + > ${OUTDIR}/${PREFIX}\_plus.bedGraph
bedtools genomecov -bg -i ${OUTDIR}/${PREFIX}.bed.gz -g ${CHINFO} -strand - > ${OUTDIR}/${PREFIX}\_minus.noinv.bedGraph

## Invert minus strand.
echo "Inverting minus strand..."
cat ${OUTDIR}/${PREFIX}\_minus.noinv.bedGraph | awk 'BEGIN{OFS="\t"} {print $1,$2,$3,-1*$4}' > ${OUTDIR}/${PREFIX}\_minus.bedGraph ## Invert read counts on the minus strand.

# Sort and convert to bigWig
echo "Sorting and converting to bigWig..."
LC_COLLATE=C sort -k1,1 -k2,2n ${OUTDIR}/${PREFIX}\_plus.bedGraph > ${OUTDIR}/${PREFIX}\_sorted_plus.bedGraph
bedGraphToBigWig ${OUTDIR}/${PREFIX}\_sorted_plus.bedGraph ${CHINFO} ${OUTDIR}/${PREFIX}_plus.bw

LC_COLLATE=C sort -k1,1 -k2,2n ${OUTDIR}/${PREFIX}\_minus.bedGraph > ${OUTDIR}/${PREFIX}\_sorted_minus.bedGraph
bedGraphToBigWig ${OUTDIR}/${PREFIX}\_sorted_minus.bedGraph ${CHINFO} ${OUTDIR}/${PREFIX}_minus.bw

# Generate log scale bedGraphs
cat ${OUTDIR}/${PREFIX}\_minus.bedGraph | awk 'BEGIN{OFS="\t"} {print $1,$2,$3,-1*log(-1*$4)/log(10)}' > ${OUTDIR}/${PREFIX}\_log10_minus.bedGraph ## Invert read counts on the minus strand and take log
cat ${OUTDIR}/${PREFIX}\_plus.bedGraph | awk 'BEGIN{OFS="\t"} {print $1,$2,$3,log($4)/log(10)}' > ${OUTDIR}/${PREFIX}\_log10_plus.bedGraph

# Generate log scale bigWigs
echo "Generating log10 scale bigWigs..."
LC_COLLATE=C sort -k1,1 -k2,2n  ${OUTDIR}/${PREFIX}\_log10_plus.bedGraph > ${OUTDIR}/${PREFIX}\_log10_sorted_plus.bedGraph
bedGraphToBigWig ${OUTDIR}/${PREFIX}\_log10_sorted_plus.bedGraph ${CHINFO} ${OUTDIR}/${PREFIX}_log10_plus.bw

LC_COLLATE=C sort -k1,1 -k2,2n ${OUTDIR}/${PREFIX}\_log10_minus.bedGraph > ${OUTDIR}/${PREFIX}\_log10_sorted_minus.bedGraph
bedGraphToBigWig ${OUTDIR}/${PREFIX}\_log10_sorted_minus.bedGraph ${CHINFO} ${OUTDIR}/${PREFIX}_log10_minus.bw

# Generate merged bigWIgs (both strands)
bigWigMerge ${OUTDIR}/${PREFIX}\_minus.bw ${OUTDIR}/${PREFIX}\_plus.bw ${OUTDIR}/${PREFIX}\_merged.bedGraph
LC_COLLATE=C sort -k1,1 -k2,2n ${OUTDIR}/${PREFIX}\_merged.bedGraph > ${OUTDIR}/${PREFIX}\_merged_sorted.bedGraph
bedGraphToBigWig ${OUTDIR}/${PREFIX}\_merged_sorted.bedGraph ${CHINFO} ${OUTDIR}/${PREFIX}\_merged.bw

# log10 scale
bigWigMerge ${OUTDIR}/${PREFIX}\_log10_minus.bw ${OUTDIR}/${PREFIX}\_plus.bw ${OUTDIR}/${PREFIX}\_merged_log10.bedGraph
LC_COLLATE=C sort -k1,1 -k2,2n ${OUTDIR}/${PREFIX}\_merged_log10.bedGraph > ${OUTDIR}/${PREFIX}\_merged_sorted_log10.bedGraph
bedGraphToBigWig ${OUTDIR}/${PREFIX}\_merged_sorted_log10.bedGraph ${CHINFO} ${OUTDIR}/${PREFIX}\_merged_log10.bw

# Remove tmp files
echo "Removing tmp files..."
rm ${OUTDIR}/${PREFIX}\_minus.noinv.bedGraph
rm ${OUTDIR}/${PREFIX}\_minus.bedGraph
rm ${OUTDIR}/${PREFIX}\_plus.bedGraph
rm ${OUTDIR}/${PREFIX}\_sorted_minus.bedGraph
rm ${OUTDIR}/${PREFIX}\_sorted_plus.bedGraph
rm ${OUTDIR}/${PREFIX}\_log10_minus.bedGraph
rm ${OUTDIR}/${PREFIX}\_log10_plus.bedGraph
rm ${OUTDIR}/${PREFIX}\_log10_sorted_minus.bedGraph
rm ${OUTDIR}/${PREFIX}\_log10_sorted_plus.bedGraph
rm ${OUTDIR}/${PREFIX}\_merged.bedGraph
rm ${OUTDIR}/${PREFIX}\_merged_sorted.bedGraph
rm ${OUTDIR}/${PREFIX}\_merged_log10.bedGraph
rm ${OUTDIR}/${PREFIX}\_merged_sorted_log10.bedGraph

echo "Done."
