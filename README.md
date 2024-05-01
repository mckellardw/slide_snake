# `slide_snake`
![slide_snake](images/slide_snake_logo.png)
## Flexible processing for spatial transcriptomics data

# [**UNDER CONSTRUCTION**]

The goal of this project is to build a snakemake workflow for assessing different preprocessing, alignment, and quantification configurations, and digging into artifacts.  


### Install w/ `mamba`/`conda` [recommended]:
```
mamba create --name slide_snake --file=envs/slide_snake.yml
mamba activate slide_snake
```
**Currently adding smaller, rule-specific conda environments to simplify builds**

I had solve issues putting snakemake into the same conda environment- if you don't have snakemake installed locally, I would recommend making a second environment which just has snakemake in it, then loading that environment before `slide_snake`. Example (assuming you already built the above `slide_snake` env):
```
mamba create --name snakemake_only -c bioconda snakemake
mamba activate snakemake_only
mamba activate slide_snake
```  


## Quick start 
- Please see `documentation/*.md` for detailed info on custom recipes, pipeline details, etc.
- How to set up a run:
  1. Install and ensure that the executable paths are functioning (see "EXEC" in `configs/config.yaml`)
  2. Build a sample sheet, containing details on each sample you would like to analyze. 
      - Fill out the file paths for the input .fastqs (can pass multiple sequencing runs), reference genomes, etc.
  3. Comment out any unwanted output files in the target rule in `Snakefile`
  4. Run the snakemake pipeline!

## Runtime details
### Example run w/ `slurm`:
```
snakemake --cluster-config config/slurm.yaml --cluster "sbatch -p {cluster.partition} -t {cluster.time} -N {cluster.nodes} --mem {cluster.mem} -o {cluster.output} --cpus-per-task={cluster.threads}" -j 16 -k -p --nt --cluster-cancel scancel --rerun-incomplete --use-conda  --conda-frontend mamba
```

### Example run w/out `slurm`:
```
snakemake -k -p -j 32 
```

## **Helpful links:**
- [Barcode download from Curio](https://curiobioscience.com/support/barcode/)
- Extract DNB barcode whitelist for StereoSeq with [ST_BarcodeMap](https://github.com/STOmics/ST_BarcodeMap) 
  - Use the "mask format change" code mentioned in the `README`
