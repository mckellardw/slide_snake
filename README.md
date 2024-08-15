# `slide_snake`
# [***UNDER CONSTRUCTION***]
![slide_snake](images/slide_snake_logo.png)
## Flexible processing for spatial transcriptomics data

The goal of this project is to build a snakemake workflow for assessing different preprocessing, alignment, and quantification configurations, and digging into artifacts.  

## Getting started
- Please see `documentation/*.md` for detailed info on custom recipes, pipeline details, etc.

### Installation:
Install base environment w/ `mamba`/`conda`:
```
mamba env create --name slide_snake --file=envs/slide_snake.yml
mamba activate slide_snake
```

### How to set up a run:
  1. Install and ensure that the executable paths are functioning (see "EXEC" in `configs/config.yaml`)
  2. Build a sample sheet, containing details on each sample you would like to analyze. Be sure to add the path to the `SAMPLE_SHEET_PATH` variable in `configs/config.yaml`.
      - Fill out the file paths for the input .fastqs (can pass multiple sequencing runs, illumina and/or ONT data), reference genomes, etc. 
      - Add the recipe(s) for each sample, separated by spaces. 
  3. Comment out any unwanted output files in the target rule in `Snakefile`
  4. Run the snakemake pipeline!

## Runtime details
*Note*, this pipeline was written in snakemake v7
### Example local run:
```
snakemake -k -p --use-conda --conda-frontend mamba -j 56
```

### Example run w/ `slurm`:
`snakemake` version 7:
```
snakemake --cluster-config config/slurm.yaml --cluster "sbatch -p {cluster.partition} -t {cluster.time} -N {cluster.nodes} --mem {cluster.mem} -o {cluster.output} --cpus-per-task={cluster.threads}" -j 16 -k -p --nt --cluster-cancel scancel --use-conda  --conda-frontend mamba
```
*Note*, rule-specific resource usage found in `config/slurm.yaml`

`snakemake` version 8:
Make sure the [slurm plugin](https://snakemake.github.io/snakemake-plugin-catalog/plugins/executor/slurm.html) is installed first!
```
snakemake -k -p --nt --use-conda --conda-frontend mamba --executor slurm --workflow-profile profiles/slurm -j 24
```

## **Helpful links:**
- [Barcode download from Curio](https://curiobioscience.com/support/barcode/)
- Extract DNB barcode whitelist for StereoSeq with [ST_BarcodeMap](https://github.com/STOmics/ST_BarcodeMap) 
  - Use the "mask format change" code mentioned in the `README`
