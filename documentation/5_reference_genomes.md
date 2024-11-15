# Reference genome/transcriptome preparation

## Option #1 - `ref_snake`
I have a separate snakemake pipeline that I use to build my references, `ref_snake`, which you can find [here](https://github.com/mckellardw/ref_snake)

## Option #2 - build them yourself

### STAR
```
STAR \
    --runThreadN ${NCORES} \
    --runMode genomeGenerate \
    --genomeDir ${OUTDIR} \
    --genomeFastaFiles ${FASTA_GENOME} \
    --sjdbGTFfile ${GENES_GTF} \
    --sjdbGTFfeatureExon exon
```

### kallisto
Using `kb-python`:
```
kb ref \
    -i transcriptome.idx \
    -g t2g.txt \
    -f1 ${cDNA_FASTA} \
    ${GENOME_FASTA} ${GENES_GTF}
```