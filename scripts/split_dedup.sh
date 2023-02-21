#!/usr/bin/bash
# split_dedup.sh - a bash script that deduplicates .bam files, aligned with STARsolo
#                - deduplicates chr-by-chr, to reduce run time and memory requirements
# Usage:
# bash chr_split_dedup.sh /path/to/input.bam /path/to/whitelist num_cores /path/to/output/directory /path/to/tmp/directory file.LOG

# bash /home/dwm269/DWM_utils/seq_utils/chr_split_dedup.sh ds.bam /workdir/dwm269/totalRNA/STRS-HD/data/STARsolo/SH2A/tmp/whitelist.txt 24 ./OUTPUT_DEDUP DEDUP.LOG

# chr_split_dedup (v1.0) - input a .bam file, output the deduplicated .bam file in OUTDIR

#TODO: python script instead of bash to make parallelization a bit easier?

INBAM=$1 # path to .bam file (sorted & indexed already!)
BB=$2 # path to barcode whitelist
CORE=$3 # number of cores for parallelization
OUTBAM=$4 # output/deduped bam path
TMPDIR=$5
LOG=$6


mkdir -p ${TMPDIR}
cd ${TMPDIR}

echo "Log for chr_split_dedup:" > ${LOG}
echo "Log info can be found in "${LOG}

PREFIX=`echo ${INBAM} | rev | cut -d / -f 1 | cut -d . -f 2- | rev`

echo "Using "${PREFIX}" as file prefix..." >> ${LOG}

# OUTBAM=${OUTDIR}/${PREFIX}_dedup.bam

echo ".bam file location:    "${INBAM} >> ${LOG}
echo "Max cores:             "${CORE} >> ${LOG}
echo "Output location:       "${OUTBAM} >> ${LOG}
echo >> ${LOG}
echo >> ${LOG}

# Check params...
if [ ! -f ${INBAM} ]
then
    echo "Can't find ${INBAM}"
fi

# Check for .bam index file
if [ ! -f ${INBAM}.bai ]
then
    echo "Indexing input .bam file because you forgot to..."
    samtools index -@ ${CORE} ${INBAM}
fi

# Remove reads that don't have a barcode (CB)
echo "Removing reads without 'CB' or 'UB' tags..." >> ${LOG}
date >> ${LOG}

samtools view -1 -b \
-@ ${CORE} \
--tag-file CB:${BB} \
-F UB:Z:- \
${INBAM} \
> ${TMPDIR}/filter.bam

echo "Done." >> ${LOG}
echo >> ${LOG}

# split bam by chromosome/strand
echo "Splitting by chromosome..." >> ${LOG}
date >> ${LOG}
bamtools split \
-in ${TMPDIR}/filter.bam \
-reference
echo "Done." >> ${LOG}
echo >> ${LOG}

BAMLIST=${TMPDIR}/BAMLIST.tmp
ls *REF_*.bam > ${BAMLIST}

# index split .bam's
echo "Indexing split .bam files..." >> ${LOG}
date >> ${LOG}
while read BAM; do
  samtools index \
  -@ ${CORE} \
  ${BAM}
done < ${BAMLIST}
echo "Done." >> ${LOG}
echo >> ${LOG}

# dedup resulting bams one-by-one
echo "Deduplicating split .bam files..." >> ${LOG}
date >> ${LOG}
#TODO: parallelize
while read BAM; do
  samtools index -@ ${CORE} ${BAM}
  umi_tools dedup \
  -I ${BAM} \
  --extract-umi-method=tag \
  --umi-tag=UB \
  --cell-tag=CB \
  --method=unique \
  --per-cell \
  --unmapped-reads=discard \
  --log ${OUTDIR}/umitools.log \
  -S ${TMPDIR}/dedup_${BAM}
done < ${BAMLIST}
# --output-stats={params.OUTPUT_PREFIX} \
echo "Done." >> ${LOG}
echo >> ${LOG}

# merge, sort, and index dedup'ed .bams
echo "Merging, sorting, and indexing deduplicated .bam files..." >> ${LOG}
date >> ${LOG}
samtools merge \
-f \
-o ${TMPDIR}/dedup_merge.bam \
-@ ${CORE} \
--write-index \
${TMPDIR}/dedup_*REF_*.bam

samtools sort ${TMPDIR}/dedup_merge.bam > ${OUTBAM}
samtools index -@ ${CORE} ${OUTBAM} && rm -r ${TMPDIR}
echo "Done." >> ${LOG}
echo >> ${LOG}
