# Installing `slide_snake`

## Install from provided .yml file:
Install base environment w/ `mamba`/`conda`:
```
mamba env create --name slide_snake --file=envs/slide_snake.yml
mamba activate slide_snake
```

## Install from scratch
Install base environment w/ `mamba`/`conda`:
```
mamba create --name slsn \
    bioconda::samtools bioconda::pysam bioconda::star \
    conda-forge::pandas conda-forge::pigz conda-forge::editdistance

mamba activate slsn

mamba install bioconda::snakemake
pip install biopython
```
