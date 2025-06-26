# Reference Genome and Transcriptome Preparation for `slide_snake`

High-quality reference genomes and transcriptomes are essential for accurate spatial RNA-seq analysis. This guide covers the preparation of all necessary reference files for `slide_snake`, including alignment indices, gene annotations, and platform-specific resources.

## Overview

`slide_snake` requires several types of reference files:
- **STAR genome indices** for genome alignment
- **Kallisto transcriptome indices** for pseudoalignment
- **Gene annotations** (GTF/GFF format) for quantification
- **Transcript sequences** (FASTA format) for ONT analysis
- **rRNA references** for contamination filtering
- **Platform-specific whitelists** for barcode correction

## Quick Start Options

### Option 1: Use `ref_snake` (Recommended)

`ref_snake` is a companion Snakemake pipeline that automates reference preparation:

```bash
# Clone ref_snake
git clone https://github.com/mckellardw/ref_snake.git
cd ref_snake

# Configure for your species
# Edit config.yaml with desired genome versions

# Run reference preparation
snakemake -j 8 --use-conda
```

**Advantages:**
- Automated download and processing
- Version control for reproducibility
- Handles multiple species simultaneously
- Generates all required file types
- Includes quality control steps

### Option 2: Manual Preparation

For custom genomes or specific requirements, build references manually following the instructions below.

## Reference File Requirements

### Core Files (Required)

| File Type | Format | Description | Used By |
|-----------|--------|-------------|---------|
| Genome FASTA | `.fa/.fasta` | Reference genome sequence | STAR, minimap2 |
| Gene Annotations | `.gtf/.gff` | Gene structure annotations | STAR, quantification |
| STAR Index | Directory | STAR alignment index | STARsolo |
| Kallisto Index | `.idx` | Kallisto pseudoalignment index | kb-python |
| Transcript-to-Gene Map | `.txt` | Transcript ID to gene ID mapping | kb-python |

### Optional Files

| File Type | Format | Description | Used By |
|-----------|--------|-------------|---------|
| Transcript FASTA | `.fa/.fasta` | Transcript sequences | ONT minimap2 |
| rRNA Reference | `.fa/.fasta` | Ribosomal RNA sequences | rRNA filtering |
| rRNA Annotations | `.gtf` | rRNA gene annotations | rRNA quantification |

## Species-Specific Configurations

### Human (Homo sapiens)

#### GENCODE (Recommended)
```bash
# Genome and annotations
GENOME_VERSION="GRCh38"
GENCODE_VERSION="47"

# URLs
GENOME_URL="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_${GENCODE_VERSION}/GRCh38.primary_assembly.genome.fa.gz"
GTF_URL="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_${GENCODE_VERSION}/gencode.v${GENCODE_VERSION}.annotation.gtf.gz"
TRANSCRIPTS_URL="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_${GENCODE_VERSION}/gencode.v${GENCODE_VERSION}.transcripts.fa.gz"
```

#### ENSEMBL (Alternative)
```bash
# Genome and annotations  
ENSEMBL_VERSION="112"
GENOME_URL="https://ftp.ensembl.org/pub/release-${ENSEMBL_VERSION}/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"
GTF_URL="https://ftp.ensembl.org/pub/release-${ENSEMBL_VERSION}/gtf/homo_sapiens/Homo_sapiens.GRCh38.${ENSEMBL_VERSION}.gtf.gz"
```

### Mouse (Mus musculus)

#### GENCODE
```bash
# Genome and annotations
GENOME_VERSION="GRCm39"
GENCODE_VERSION="M36"

# URLs
GENOME_URL="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_${GENCODE_VERSION}/GRCm39.primary_assembly.genome.fa.gz"
GTF_URL="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_${GENCODE_VERSION}/gencode.v${GENCODE_VERSION}.annotation.gtf.gz"
TRANSCRIPTS_URL="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_${GENCODE_VERSION}/gencode.v${GENCODE_VERSION}.transcripts.fa.gz"
```

### Other Species

For other species, use ENSEMBL:
```bash
# Replace 'species_name' with your species (e.g., rattus_norvegicus, danio_rerio)
SPECIES="species_name"
ENSEMBL_VERSION="112"

GENOME_URL="https://ftp.ensembl.org/pub/release-${ENSEMBL_VERSION}/fasta/${SPECIES}/dna/"
GTF_URL="https://ftp.ensembl.org/pub/release-${ENSEMBL_VERSION}/gtf/${SPECIES}/"
```

## Building Reference Indices

### STAR Genome Index

```bash
# Set paths
GENOME_FA="/path/to/genome.fa"
GENES_GTF="/path/to/genes.gtf"
STAR_INDEX_DIR="/path/to/STAR_index"
THREADS=16

# Create index directory
mkdir -p ${STAR_INDEX_DIR}

# Build STAR index
STAR \
    --runThreadN ${THREADS} \
    --runMode genomeGenerate \
    --genomeDir ${STAR_INDEX_DIR} \
    --genomeFastaFiles ${GENOME_FA} \
    --sjdbGTFfile ${GENES_GTF} \
    --sjdbGTFfeatureExon exon \
    --sjdbOverhang 100
```

**Important Parameters:**
- `--sjdbOverhang`: Set to (read_length - 1), or use 100 for variable lengths
- `--genomeSAindexNbases`: May need adjustment for small genomes

### Kallisto Transcriptome Index

Using `kb-python` (recommended):
```bash
# Set paths
GENOME_FA="/path/to/genome.fa"
GENES_GTF="/path/to/genes.gtf"
KB_INDEX="/path/to/transcriptome.idx"
T2G_FILE="/path/to/t2g.txt"

# Build index with kb ref
kb ref \
    -i ${KB_INDEX} \
    -g ${T2G_FILE} \
    -f1 transcripts.fa \
    ${GENOME_FA} ${GENES_GTF}
```

Using `kallisto` directly:
```bash
# Extract transcriptome
gffread -w transcripts.fa -g ${GENOME_FA} ${GENES_GTF}

# Build kallisto index
kallisto index \
    -i ${KB_INDEX} \
    transcripts.fa

# Create transcript-to-gene mapping
# (Custom script required)
```

### Minimap2 Indices (ONT)

For genome alignment:
```bash
# Build minimap2 index
minimap2 -d genome.mmi ${GENOME_FA}
```

For transcriptome alignment:
```bash
# Build transcriptome index
minimap2 -d transcriptome.mmi ${TRANSCRIPTS_FA}
```

## rRNA Reference Preparation

### Download rRNA Sequences

```bash
# Human rRNA (example)
# 18S rRNA
wget "https://www.ncbi.nlm.nih.gov/sviewer/viewer.fcgi?id=555853&db=nuccore&report=fasta" -O 18S.fa

# 28S rRNA  
wget "https://www.ncbi.nlm.nih.gov/sviewer/viewer.fcgi?id=555858&db=nuccore&report=fasta" -O 28S.fa

# 5.8S rRNA
wget "https://www.ncbi.nlm.nih.gov/sviewer/viewer.fcgi?id=555857&db=nuccore&report=fasta" -O 5.8S.fa

# 5S rRNA
wget "https://www.ncbi.nlm.nih.gov/sviewer/viewer.fcgi?id=555856&db=nuccore&report=fasta" -O 5S.fa

# Combine all rRNA sequences
cat 18S.fa 28S.fa 5.8S.fa 5S.fa > rRNA_combined.fa
```

### Build rRNA Index

For BWA-based filtering:
```bash
# Build BWA index
bwa index rRNA_combined.fa
```

## Platform-Specific Resources

### Barcode Whitelists

#### 10x Genomics Visium
```bash
# Download from 10x website
wget https://cf.10xgenomics.com/supp/spatial-exp/visium-v1-spatial-barcodes.txt

# Or use provided whitelist in slide_snake
# resources/visium_whitelist/whitelist.txt
```

#### Curio Seeker
```bash
# Download from Curio website
# https://curiobioscience.com/support/barcode/
# Or use custom whitelist based on your array
```

#### DBIT-seq
```bash
# Use provided whitelists in slide_snake
# resources/dbit_whitelist/Spatial_barcode_50x50.txt
# resources/dbit_whitelist/Spatial_barcode_100x100.txt
```

#### StereoSeq/STOmics
```bash
# Use ST_BarcodeMap to extract from chip file
git clone https://github.com/STOmics/ST_BarcodeMap.git
cd ST_BarcodeMap

# Extract barcode whitelist
ST_BarcodeMap-0.0.1 \
    --in chip_file.h5 \
    --out whitelist.txt \
    --action 3 \
    --thread 8
```

## Quality Control and Validation

### Validate STAR Index
```bash
# Check STAR index
STAR \
    --genomeDir ${STAR_INDEX_DIR} \
    --genomeLoad LoadAndExit

# Should show genome loading without errors
```

### Validate Kallisto Index
```bash
# Check kallisto index
kallisto inspect ${KB_INDEX}

# Should show transcriptome statistics
```

### Test Alignments
```bash
# Test with a small FASTQ file
STAR \
    --genomeDir ${STAR_INDEX_DIR} \
    --readFilesIn test_R1.fq.gz test_R2.fq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --outFileNamePrefix test_
```

## Storage and Organization

### Recommended Directory Structure
```
/path/to/references/
├── homo_sapiens/
│   ├── GRCh38_GENCODE47/
│   │   ├── genome/
│   │   │   ├── genome.fa
│   │   │   ├── genes.gtf
│   │   │   └── transcripts.fa
│   │   ├── STAR/
│   │   │   └── [STAR index files]
│   │   ├── kallisto/
│   │   │   ├── transcriptome.idx
│   │   │   └── t2g.txt
│   │   └── rRNA/
│   │       ├── rRNA.fa
│   │       └── rRNA.gtf
├── mus_musculus/
│   └── GRCm39_GENCODEM36/
│       └── [similar structure]
└── platform_resources/
    ├── visium_whitelist/
    ├── dbit_whitelist/
    └── stomics_whitelist/
```

### Version Control
- Use consistent naming with version numbers
- Document source and date of download
- Keep checksums for validation
- Use symbolic links for "current" versions

## Sample Sheet Configuration

Once references are prepared, configure your sample sheet:

```text
sampleID,species,STAR_ref,genes_gtf,kb_idx,kb_t2g,genome_fa,cDNA_fa,rRNA_ref,rRNA_gtf
sample1,human,/ref/homo_sapiens/STAR/,/ref/homo_sapiens/genes.gtf,/ref/homo_sapiens/kallisto/transcriptome.idx,/ref/homo_sapiens/kallisto/t2g.txt,/ref/homo_sapiens/genome.fa,/ref/homo_sapiens/transcripts.fa,/ref/homo_sapiens/rRNA/rRNA.fa,/ref/homo_sapiens/rRNA/rRNA.gtf
```

## Troubleshooting

### Common Issues

**STAR index build fails:**
- Check genome file format (FASTA headers should be simple)
- Ensure sufficient memory (32GB+ for human genome)
- Verify GTF format compatibility

**Kallisto index problems:**
- Ensure transcript sequences match GTF annotations
- Check for duplicate transcript IDs
- Validate FASTA format

**Memory issues:**
- Use `--genomeSAindexNbases` parameter for STAR
- Consider splitting large genomes
- Use SSD storage for faster I/O

### Performance Tips

**Faster index building:**
- Use SSD storage
- Maximize available RAM
- Use multiple threads
- Pre-sort input files

**Storage optimization:**
- Compress unused files
- Use hard links for duplicate files
- Regular cleanup of temporary files

## Best Practices

1. **Consistency**: Use the same genome version across all analyses
2. **Documentation**: Keep detailed records of sources and versions
3. **Validation**: Always test indices with small datasets first
4. **Backup**: Maintain backups of reference files
5. **Updates**: Regularly update to newer genome versions
6. **Sharing**: Use shared storage for team projects

This comprehensive guide ensures you have high-quality, properly formatted reference files for accurate spatial RNA-seq analysis with `slide_snake`.