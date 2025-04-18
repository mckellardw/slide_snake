# Toy data for Patho-dBIT

## `Human_Lymphoma_DLBCL_Region2_10um` 
- [GEO link](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM8454085)
- [github](https://github.com/Zhiliang-Bai/Patho-DBiT)
- [barcode map](https://github.com/mckellardw/slide_snake/blob/main/resources/dbit_whitelist/Spatial_barcode_100x100.txt)
- [Citation](https://pubmed.ncbi.nlm.nih.gov/39353436/)
- Patho-dBIT data generated on Diffuse Large B-Cell Lymphoma 
  -   10um pixels
  -   100-by-100 in situ barcoding schema
- Only Illumina reads available in toy data (10,000 reads)

How to download:
- Update the file paths in the bash script below to your local machine
- Run the bash script [`download.sh`] to download the data
```bash
DATADIR="/path/to/data/patho_dbit"
mkdir -p ${DATADIR}
cd ${DATADIR}

ACC_LIST="/path/to/slide_snake/data/test_pathodbit/SRR_Acc_List.txt"

NTHREADS=24

## load in SRR IDs
readarray -t SRR < ${ACC_LIST}

for i in "${SRR[@]}"
do
  echo ${i}
  fastq-dl \
    --accession ${i} \
    --outdir ${DATADIR} \
    --prefix ${i} \
    --cpus ${NTHREADS} \
    --group-by-sample
done
```
