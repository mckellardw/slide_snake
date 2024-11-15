# Installing `slide_snake`

## Install from provided .yml file:
Install base environment w/ `mamba`/`conda`:
```
mamba env create --file=envs/slsn.yml
mamba activate slsn
```

## Install from scratch
Install base environment w/ `mamba`/`conda`:
```
mamba create --name slsn \
    bioconda::snakemake bioconda::samtools bioconda::pysam  \
    conda-forge::pandas conda-forge::pigz 

mamba activate slsn
```
