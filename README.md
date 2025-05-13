# `slide_snake`
# [***UNDER CONSTRUCTION***]
![slide_snake](images/slide_snake_logo.png)
## Flexible processing for spatial RNA-sequencing data

`slide_snake` is a snakemake workflow for preprocessing, alignment, and quantification of spatial RNA-sequencing data, designed to work for any commercial or custom platform. 

## Getting started
- Please see `documentation/*.md` for detailed info on custom recipes, pipeline details, etc.

### Installation:
Install base environment w/ `mamba`:
```
mamba env create --name slsn --file=envs/slsn.yml
mamba activate slsn
```

### How to set up a run:
  1. Install and ensure that the executable paths are functioning
  2. Build a sample sheet, containing details on each sample you would like to analyze. Be sure to add the path to the `SAMPLE_SHEET_PATH` variable in `configs/config.yaml`.
      - Fill out the file paths for the input .fastqs (can pass multiple sequencing runs, illumina and/or ONT data), reference genomes, etc. 
      - Add the recipe(s) for each sample, separated by spaces. 
  3. Comment out any unwanted output files in the target rule in `Snakefile`
  4. Run the snakemake pipeline!

## Runtime details
### Example local run:
```
snakemake -k -p --use-conda --conda-frontend mamba -j 56
```

### Example run w/ `slurm`:
Make sure the [slurm plugin](https://snakemake.github.io/snakemake-plugin-catalog/plugins/executor/slurm.html) is installed first! (it is included in `envs/slsn.yml`)
```
snakemake -k -p --nt --use-conda --conda-frontend mamba --executor slurm --workflow-profile profiles/slurm -j 24
```
*Note*, this pipeline was written with snakemake v8

Or, my preferred way to run it in the background:
```
nohup snakemake -k -p --use-conda --conda-frontend mamba -j 36 > snake.log 2>&1 &
```


# Helpful links:
- [Barcode download from Curio](https://curiobioscience.com/support/barcode/)
- Extract DNB barcode whitelist for StereoSeq with [ST_BarcodeMap](https://github.com/STOmics/ST_BarcodeMap) 
  - Use the "mask format change" code mentioned in the `README`; example:
  ```
  ST_BarcodeMap-0.0.1 --in B01807A3.barcodeToPos.h5 --out B01807A3.barcodeToPos.txt --action 3 --thread 24
  ```
- "Documentation" on `kallisto-lr` [link](https://github.com/pachterlab/kallisto/issues/456)