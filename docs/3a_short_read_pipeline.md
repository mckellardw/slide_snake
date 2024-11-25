# Short read preprocessing info

## Trimming

### Five Prime Adapter(s)
- Reverse complement of the SlideSeq TSO to reduce strand invasion artifacts: `CCCATTCACTCTGCGTTGATACCAGCTT`

### Three Prime Adapter(s)
- "A" homopolymers: `100-A`
- "G" homopolymers, important for Illumina 2-color sequencers: `100-G`
- "T" homopolymers: `100-T`
- Nextera adapter sequence: `CTGTCTCTTATA`
- Reverse complement of Nextera sequence: `TATAAGAGACAG`
- Curio template switch oligo (TSO) - remove any polyadenylated TSOs: `AAGCTGGTATCAACGCAGAGTGAATGGG`
- Curio R1 internal adapter - shows up in R2, and used for R1 trimming: `TCTTCAGCGTTCCCGAGA`
- Reverse of Curio R1 adapter: `AGAGCCCTTGCGACTTCT`
- Illumina universal sequence: `AGATCGGAAGAG`

## Alignment

### STAR
- After adapter trimming, reads are aligned with `STARsolo` and `kallisto`/`bustools` to generate count matrices
- Outputs are in `{SAMPLE_ID}/STARsolo/Solo.out` & `{SAMPLE_ID}/kb/counts_unfiltered`
- Different recipes are written out in `resources/recipe_sheet.csv`, and must be specified for each sample within the sample sheet

#### Generating References

##### STAR Reference
This is a typical STAR reference that you would use for any other alignment job. Here is an example code snippet:

```bash
FASTA_DIR="/path/to/GENCODE_M31/GRCm39.genome.fa"
GENES_DIR="/path/to//GENCODE_M31/gencode.vM31.annotation.gtf"

OUTDIR="/workdir/dwm269/genomes/mm39_all/STAR_GRCm39_GENCODEM31"

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

## Barcode Handling

### STAR
- Removed the linker sequence in R1 so that the `1MM_multi` barcode correction in `STARsolo` can be used
- Barcode & UMI parameters for [`STAR`](https://github.com/alexdobin/STAR/blob/master/docs/STARsolo.md):
```bash
--soloType CB_UMI_Simple \
--soloUMIstart 14 \
--soloUMIlen 7 \
--soloCBstart 1 \
--soloCBlen 14
```

### kallisto/bustools
- Barcode & UMI parameters for `kallisto bus`:
```bash
-x 0,0,14:0,14,21:1,0,0
```

## QC
- `fastqc` is run on R2 files before (`preTrim`) and after (`postTrim`) adapter trimming
- `qualimap` is used to assess `STAR` alignment

## Unmapped Read Analysis
- `fastqc` is run on all unmapped reads from `STAR`
- The most abundant unmapped reads are also aligned with `blast`

## Small/Micro RNA Analysis
**Work-in-progress...**
- Currently have rules set up to align to databases listed below (with `bowtie2`) and generate a count matrix (note, this uses `STAR`'s cell/bead/spot barcoding)
- Need to better optimize alignment params

### miRge3.0
- [Link to documentation](https://mirge3.readthedocs.io/en/latest/quick_start.html)
- [Link to library download](https://sourceforge.net/projects/mirge3/files/miRge3_Lib/)
- Bulk miRg3.0 is implemented, but it is quite slow. I recommend commenting out the target rule(s) unless you have small/total RNA data for which you need this analysis.

### Small RNA Reference Databases
- **piRNAs** - [pirbase](http://bigdata.ibp.ac.cn/piRBase/), v3.0
  - *Note* - use "gold standard" piRs because there are way too many predicted sequences...
- **miRNAs** - [mirbase](https://mirbase.org/), v22
  - Use both hairpin & mature sequences, in that order