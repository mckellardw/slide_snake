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

## Reference Genome Info

### STAR Reference
This is a typical STAR reference that you would use for any other alignment job. Here is an example code snippet:

```bash
FASTA_DIR="/path/to/GENCODE_M36/GRCm39.genome.fa"
GENES_DIR="/path/to/GENCODE_M36/gencode.vM36.annotation.gtf"

OUTDIR="/workdir/dwm269/genomes/mm39_all/STAR_GRCm39_GENCODEM36"

mkdir -p ${OUTDIR}
cd ${OUTDIR}

STAR \
  --runThreadN 16 \
  --runMode genomeGenerate \
  --genomeDir ${OUTDIR} \
  --genomeFastaFiles ${FASTA_DIR} \
  --sjdbGTFfile ${GENES_DIR} \
  --sjdbGTFfeatureExon exon
```

*You can find the reference files on [GENCODE's website](https://www.gencodegenes.org/mouse/)*

### kallisto/BUStools
TODO
