# Test Datasets

This directory contains toy datasets for testing various spatial transcriptomics technologies supported by slide_snake. Each dataset is reduced in size for quick testing and demonstration purposes.

## Visium/STRS (`test_visium/`)

### STRS/Skeletal muscle - "vy3C"
- [GEO link](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE161318)
- [barcode map](https://github.com/mckellardw/slide_snake/blob/main/resources/visium_whitelist/visium-v1_coordinates.txt)
- [Citation](https://www.nature.com/articles/s41587-022-01517-6) - see Fig 2
- Visium/STRS data
- Both Illumina and ONT reads available in toy data

### Visium/olfactory bulb - "sit_mob"
- [GEO link](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE153859)
- [barcode map](https://github.com/mckellardw/slide_snake/blob/main/resources/visium_whitelist/visium-v1_coordinates.txt)
- [Citation](https://doi.org/10.1093/nar/gkad169)
- First Visium long-read (Oxford Nanopore) data
- Both Illumina and ONT reads available in toy data (#TODO)


## Seeker/SlideSeq_v2 (`test_seeker/`)

### `Heart_Control` - Seeker/heart
- [GEO link](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE161318)
- [barcode map](https://github.com/mckellardw/slide_snake/blob/main/data/test_seeker/A0004_043_BeadBarcodes.txt)
- [Citation](https://www.nature.com/articles/s44161-022-00138-1) - see Fig 4
- Seeker data generated on reovirus-infected mouse heart
- Both Illumina and ONT reads available in toy data

### Seeker/ONT #TODO

- [NCBI link](https://dataview.ncbi.nlm.nih.gov/object/PRJNA1116561?reviewer=uc7v6e10toqp4lllatmu6vvf53) - still not public, will add once they do...

## Patho-dBIT (`test_pathodbit/`)

### `Human_Lymphoma_DLBCL_Region2_10um` 
- [GEO link](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM8454085)
- [github](https://github.com/Zhiliang-Bai/Patho-DBiT)
- [barcode map](https://github.com/mckellardw/slide_snake/blob/main/resources/dbit_whitelist/Spatial_barcode_100x100.txt)
- [Citation](https://pubmed.ncbi.nlm.nih.gov/39353436/)
- Patho-dBIT data generated on Diffuse Large B-Cell Lymphoma 
  - 10um pixels
  - 100-by-100 in situ barcoding schema
- Only Illumina reads available in toy data (10,000 reads)

### How to download Patho-dBIT data:
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

## DecoderSeq (`test_decoderseq/`)

- [GEO link](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE2358968)
- [barcode map](https://github.com/mckellardw/slide_snake/blob/main/resources/decoderseq_whitelist/barcodeslist.txt)
- [Citation](https://www.nature.com/articles/s41587-023-02086-y) - see Fig 2
- DecoderSeq data

## MAGIC-seq (`test_magicseq/`)

### `TODO` 
- [GEO link](TODO)
- [github](https://github.com/bioinfo-biols/MAGIC-seq)
- [doi](https://doi.org/10.1038/s41588-024-01906-4) 
- [barcode map](TODO)
- [Citation](TODO)
- SAMPLE DESCRIPTION TODO
- Only BGI/MGI/Complete Genomics reads available in toy data (10,000 reads)

## microST (`test_microST/`)

*TODO - Dataset information to be added*

---

## Sample Sheets

The data directory also contains example sample sheets for testing:
- `test_sample_sheet_minimal.csv` - Minimal sample sheet format
- `test_sample_sheet.csv` - Complete sample sheet with all optional fields
