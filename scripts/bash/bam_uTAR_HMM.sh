# bash scripts/bash/bam_uTAR_HMM.sh --bam INPUT_BAM --threads N_THREADS --mem MEM --mergebp MERGEBP --thresh THRESH --outdir OUTDIR

# Default values
N_THREADS=5
MERGEBP=500
THRESH=10000000
MINCOV=5

# Parse arguments
while [[ "$#" -gt 0 ]]; do
    case $1 in
        --bam) INPUT_BAM="$2"; shift ;;
        --threads) N_THREADS="$2"; shift ;;
        --mem) MEM="$2"; shift ;;
        --mergebp) MERGEBP="$2"; shift ;;
        --thresh) THRESH="$2"; shift ;;
        --outdir) OUTDIR="$2"; shift ;;
        *) echo "Unknown parameter: $1"; exit 1 ;;
    esac
    shift
done

# Validate required arguments
if [[ -z "$INPUT_BAM" || -z "$OUTDIR" ]]; then
    echo "Error: --bam and --outdir are required arguments."
    exit 1
fi

# Set up directories and variables
PREFIX="HMMout"
TMPDIR=${OUTDIR}/${PREFIX}
mkdir -p ${TMPDIR}

# Ensure the output directory exists
mkdir -p "$OUTDIR"

# Function to log messages with timestamps
log_with_timestamp() {
    echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"
}

# Log input parameters
echo "input .BAM                  ${INPUT_BAM}"
echo "output directory:           ${OUTDIR}"
echo "tmp directory:              ${TMPDIR}"
echo "number of threads:          ${N_THREADS}"
echo "memory usage:               ${MEM}"
echo "minimum coverage:           ${MINCOV}"
echo "thresholded at 1 in ${THRESH} reads"
echo ""
echo "Reads spanning over splicing junction will join HMM blocks"
echo "To avoid that, split reads into small blocks before input to groHMM"
echo "Splitting and sorting reads..."

bedtools bamtobed -i ${INPUT_BAM} -split \
| LC_ALL=C sort -k1,1V -k2,2n --buffer-size=${MEM} --parallel=${N_THREADS} \
| gzip > ${TMPDIR}/${PREFIX}_split.sorted.bed.gz

log_with_timestamp "Splitting and sorting completed."

zcat ${TMPDIR}/${PREFIX}_split.sorted.bed.gz | awk -v TMPDIR=${TMPDIR} '{print $0 >> TMPDIR"/"$1".bed"}'
find ${TMPDIR} -name "*.bed" -size -1024k -delete
wc ${TMPDIR}/chr*.bed -l > ${TMPDIR}/chr_read_count.txt

log_with_timestamp "Running groHMM on each individual chromosome..."

wait_a_second() {
	joblist=($(jobs -p))
    while (( ${#joblist[*]} >= ${N_THREADS} ))
	    do
	    sleep 1
	    joblist=($(jobs -p))
	done
}

for f in ${TMPDIR}/*.bed
do
  wait_a_second
  log_with_timestamp "Processing chromosome file: ${f}"
  Rscript --vanilla scripts/R/uTAR_HMM.R --input_file=${f} --output_dir=${TMPDIR} 2>&1 & pids+=($!)
done
wait "${pids[@]}"

log_with_timestamp "HMM processing completed for all chromosomes."

log_with_timestamp "Merging HMM blocks within ${MERGEBP}bp..."
for f in ${TMPDIR}/*_HMM.bed
do
  LC_ALL=C sort -k1,1V -k2,2n --parallel=${N_THREADS} ${f} > ${f}.sorted.bed
  cat ${f}.sorted.bed | grep + > ${f}_plus
  cat ${f}.sorted.bed | grep - > ${f}_minus
  bedtools merge -s -d ${MERGEBP} -i ${f}_plus > ${f}_plus_merge${MERGEBP} & pids2+=($!)
  bedtools merge -s -d ${MERGEBP} -i ${f}_minus > ${f}_minus_merge${MERGEBP} & pids2+=($!)
  wait_a_second
done
wait "${pids2[@]}"

log_with_timestamp "Combining HMM output from all chromosomes..."
cat ${TMPDIR}/*_HMM.bed_plus_merge${MERGEBP}  \
| awk 'BEGIN{OFS="\t"} {print $0, ".", ".", "+"}' \
> ${TMPDIR}/${PREFIX}_merge${MERGEBP}

cat ${TMPDIR}/*_HMM.bed_minus_merge${MERGEBP} \
| awk 'BEGIN{OFS="\t"} {print $0, ".", ".", "-"}' \
>> ${TMPDIR}/${PREFIX}_merge${MERGEBP}

mkdir -p ${TMPDIR}/toremove
for f in ${TMPDIR}/*_HMM.bed
do
	mv ${f}.sorted.bed ${f}_plus ${f}_minus ${f}_plus_merge${MERGEBP} ${f}_minus_merge${MERGEBP} ${TMPDIR}/toremove/.
done

log_with_timestamp "Sorting combined .bed file..."
f=${TMPDIR}/${PREFIX}
LC_ALL=C sort -k1,1V -k2,2n ${f}_merge${MERGEBP} --parallel=${N_THREADS} > ${f}_merge${MERGEBP}.sorted.bed

zcat ${TMPDIR}/${PREFIX}_split.sorted.bed.gz | awk {'print $1'} | sort -k1,1V | uniq | sed s/$/'\t42'/ > ${TMPDIR}/tmp.genome

log_with_timestamp "Calculating the coverage..."
bedtools coverage -nonamecheck -a ${f}_merge${MERGEBP}.sorted.bed -b <(zcat ${TMPDIR}/${PREFIX}_split.sorted.bed.gz) -s -counts -split -sorted -g ${TMPDIR}/tmp.genome > ${f}_merge${MERGEBP}.sorted.bed_count
rm ${TMPDIR}/tmp.genome

log_with_timestamp "Filtering the HMM blocks by coverage..."
cat ${f}_merge${MERGEBP}.sorted.bed_count | awk 'BEGIN{OFS="\t"} ($7 >= '$MINCOV'){print $1, $2, $3, $4, $5, $6, $7}' | gzip > ${OUTDIR}/TAR_raw.bed.gz

log_with_timestamp "#### Please examine if major chromosomes are all present in the final TAR_raw.bed.gz file ####"
zcat ${OUTDIR}/TAR_raw.bed.gz | cut -f 1 | uniq

# log_with_timestamp "Linking the final TAR_reads.bed.gz file to the working directory"
# ln -sf ${TMPDIR}/TAR_reads.bed.gz ${OUTDIR}/TAR_reads.bed.gz

log_with_timestamp "Moving intermediate files to ${TMPDIR}/toremove ..."
mv ${TMPDIR}/chr* ${TMPDIR}/toremove/.

# Only gzip the final output .bed file
pigz -p ${N_THREADS} ${TMPDIR}/TAR_reads.bed &

log_with_timestamp "Done!"
