# Sample Sheets in `slide_snake`

Sample sheets are the central configuration files that tell `slide_snake` how to process your spatial RNA-seq data. They provide flexibility in specifying different analysis parameters, reference genomes, and processing recipes for each sample.

## Overview

The sample sheet is a CSV file that contains all the necessary information about your samples, including:
- Input FASTQ file paths
- Reference genome locations  
- Processing recipes to apply
- Platform-specific parameters
- Quality control settings

## Quick Start

1. **Use the template**: Start with `docs/example_sample_sheet.csv`
2. **Edit for your data**: Modify paths and parameters for your samples
3. **Update config**: Set `SAMPLE_SHEET_PATH` in `config.yaml` to point to your sample sheet
4. **Validate**: Run `snakemake -n` to check for errors

## Sample Sheet Columns

### Required Columns

| Column | Description | Example |
|--------|-------------|---------|
| `sampleID` | Unique identifier for each sample | `sample_001_visium` |
| `fastq_R1` | Path to Read 1 FASTQ file(s) | `data/sample1_R1.fq.gz` |
| `fastq_R2` | Path to Read 2 FASTQ file(s) | `data/sample1_R2.fq.gz` |
| `recipe` | Processing recipe(s) for short reads | `visium` or `visium seeker` |
| `species` | Species name for reference selection | `human`, `mouse`, `rat` |

### Platform Detection Columns

| Column | Description | Example |
|--------|-------------|---------|
| `whitelist` | Path to barcode whitelist file | `resources/visium_whitelist/whitelist.txt` |
| `recipe_ONT` | Recipe(s) for Oxford Nanopore data | `visium_ont` |

### Reference Genome Columns

| Column | Description | Required For | Example |
|--------|-------------|--------------|---------|
| `STAR_ref` | STAR genome index directory | STAR alignment | `/path/to/STAR_index/` |
| `kb_idx` | Kallisto index file | Kallisto alignment | `/path/to/transcriptome.idx` |
| `kb_t2g` | Transcript-to-gene mapping | Kallisto alignment | `/path/to/t2g.txt` |
| `genome_fa` | Reference genome FASTA | Various analyses | `/path/to/genome.fa` |
| `genes_gtf` | Gene annotation GTF file | Gene quantification | `/path/to/genes.gtf` |
| `cDNA_fa` | Transcript sequences FASTA | ONT transcriptome alignment | `/path/to/transcripts.fa` |

### Quality Control Columns

| Column | Description | Example |
|--------|-------------|---------|
| `rRNA_ref` | rRNA reference for filtering | `/path/to/rRNA_ref.fa` |
| `rRNA_gtf` | rRNA annotation GTF file | `/path/to/rRNA.gtf` |

## File Path Specifications

### Multiple Files
Specify multiple FASTQ files in several ways:
```csv
# Space-separated files
fastq_R1,fastq_R2
"file1_R1.fq.gz file2_R1.fq.gz","file1_R2.fq.gz file2_R2.fq.gz"

# Glob patterns (use with caution)
fastq_R1,fastq_R2
"data/sample*_R1.fq.gz","data/sample*_R2.fq.gz"
```

### Path Types
- **Absolute paths**: `/full/path/to/file.fq.gz`
- **Relative paths**: `data/file.fq.gz` (relative to workflow directory)
- **Compressed files**: `.gz`, `.bz2` files are automatically handled

## Recipe Specifications

### Single Recipe
```csv
recipe
visium
```

### Multiple Recipes
Use spaces to separate multiple recipes for benchmarking:
```csv
recipe
"visium visium_total"
```

This will run both standard and total RNA alignment strategies.

## Species and Reference Configuration

Common species configurations:

### Human (GRCh38)
```csv
species,STAR_ref,genes_gtf,kb_idx,kb_t2g
human,/path/to/GRCh38_STAR/,/path/to/gencode.v47.gtf,/path/to/human.idx,/path/to/human_t2g.txt
```

### Mouse (GRCm39) 
```csv
species,STAR_ref,genes_gtf,kb_idx,kb_t2g  
mouse,/path/to/GRCm39_STAR/,/path/to/gencode.vM36.gtf,/path/to/mouse.idx,/path/to/mouse_t2g.txt
```

## Platform-Specific Examples

### 10x Visium
```csv
sampleID,fastq_R1,fastq_R2,recipe,whitelist,species,STAR_ref,genes_gtf,kb_idx,kb_t2g
visium_sample,data/visium_R1.fq.gz,data/visium_R2.fq.gz,visium,resources/visium_whitelist/whitelist.txt,human,/ref/GRCh38_STAR,/ref/gencode.v47.gtf,/ref/human.idx,/ref/human_t2g.txt
```

### Curio Seeker
```csv
sampleID,fastq_R1,fastq_R2,recipe,whitelist,species,STAR_ref,genes_gtf,kb_idx,kb_t2g
seeker_sample,data/seeker_R1.fq.gz,data/seeker_R2.fq.gz,seeker_MatchLinker,data/barcodes.txt,mouse,/ref/GRCm39_STAR,/ref/gencode.vM36.gtf,/ref/mouse.idx,/ref/mouse_t2g.txt
```

### DBIT-seq
```csv
sampleID,fastq_R1,fastq_R2,recipe,whitelist,species,STAR_ref,genes_gtf,kb_idx,kb_t2g
dbit_sample,data/dbit_R1.fq.gz,data/dbit_R2.fq.gz,dbit-pretrim,resources/dbit_whitelist/Spatial_barcode_100x100.txt,human,/ref/GRCh38_STAR,/ref/gencode.v47.gtf,/ref/human.idx,/ref/human_t2g.txt
```

### Oxford Nanopore (ONT)
```csv
sampleID,fastq_R1,fastq_R2,recipe,recipe_ONT,whitelist,species,cDNA_fa,genes_gtf
ont_sample,data/ont_reads.fq.gz,,visium,visium_ont,resources/visium_whitelist/whitelist.txt,human,/ref/transcripts.fa,/ref/gencode.v47.gtf
```

## Advanced Configuration

### Custom Recipes
Create custom recipes by modifying `resources/recipe_sheet.csv`. See [Recipes Documentation](2_recipes.md) for details.

### rRNA Filtering
Enable ribosomal RNA filtering:
```csv
rRNA_ref,rRNA_gtf
/path/to/rRNA_BWA_index,/path/to/rRNA_annotations.gtf
```

### Quality Control Parameters
The pipeline automatically handles:
- FastQC analysis on raw and processed reads
- Alignment quality assessment with Qualimap
- Read mapping statistics
- Barcode correction statistics

## Validation and Troubleshooting

### Check Sample Sheet Format
```bash
# Dry run to validate sample sheet
snakemake -n --use-conda

# Check specific sample
snakemake -n out/{SAMPLE_ID}/illumina/STAR/{RECIPE}/counts_filtered/
```

### Common Issues

**Issue**: File paths not found
- **Solution**: Use absolute paths or verify relative paths from workflow directory

**Issue**: Recipe not recognized  
- **Solution**: Check recipe name exists in `resources/recipe_sheet.csv`

**Issue**: Reference files missing
- **Solution**: Verify all reference paths exist and are readable

**Issue**: Multiple files not processed
- **Solution**: Ensure proper quoting of space-separated file lists

### Best Practices

1. **Use absolute paths** when possible to avoid path resolution issues
2. **Test with small datasets** before running full analyses  
3. **Keep consistent naming** across samples for easier downstream analysis
4. **Document your recipes** if creating custom analysis parameters
5. **Validate references** ensure all reference files exist and are properly formatted
6. **Use version control** to track changes to your sample sheets

## Example Sample Sheets

### Complete Example
```csv
sampleID,fastq_R1,fastq_R2,recipe,recipe_ONT,whitelist,species,STAR_ref,genes_gtf,kb_idx,kb_t2g,rRNA_ref,rRNA_gtf,genome_fa,cDNA_fa
test_visium,data/visium_R1.fq.gz,data/visium_R2.fq.gz,visium,,resources/visium_whitelist/whitelist.txt,human,/ref/GRCh38_STAR,/ref/gencode.v47.gtf,/ref/human.idx,/ref/human_t2g.txt,/ref/rRNA_BWA,/ref/rRNA.gtf,/ref/GRCh38.fa,/ref/transcripts.fa
test_seeker,data/seeker_R1.fq.gz,data/seeker_R2.fq.gz,"seeker_MatchLinker seeker_std",,data/seeker_barcodes.txt,mouse,/ref/GRCm39_STAR,/ref/gencode.vM36.gtf,/ref/mouse.idx,/ref/mouse_t2g.txt,/ref/rRNA_BWA,/ref/rRNA.gtf,/ref/GRCm39.fa,/ref/mouse_transcripts.fa
```

This comprehensive sample sheet demonstrates multi-recipe analysis and complete reference specification.