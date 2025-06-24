# Installing `slide_snake`

`slide_snake` is a flexible Snakemake workflow for processing spatial RNA-sequencing data from various commercial and custom platforms. This guide will walk you through the installation process and initial setup.

## Prerequisites

Before installing `slide_snake`, ensure you have the following:
- A Unix-like operating system (Linux or macOS)
- [Mamba](https://mamba.readthedocs.io/) or [Conda](https://docs.conda.io/en/latest/) package manager
- Git for cloning the repository
- At least 8GB of available RAM (32GB+ recommended for large datasets)
- 50GB+ of free disk space for references and temporary files

## Quick Installation

### Step 1: Clone the Repository
```bash
git clone https://github.com/mckellardw/slide_snake.git
cd slide_snake
```

### Step 2: Create Environment from Provided File
Install the base environment using the provided environment file:
```bash
mamba env create --file=envs/slsn.yml
mamba activate slsn
```

This will install all necessary dependencies including Snakemake, alignment tools, and analysis packages.

### Step 3: Test Installation
Verify your installation by running a dry-run on the test data:
```bash
snakemake -n --use-conda
```

## Manual Installation (Alternative)

If you prefer to create the environment manually:
```bash
mamba create --name slsn \
    bioconda::snakemake \
    bioconda::samtools \
    bioconda::pysam \
    conda-forge::pandas \
    conda-forge::pigz \
    bioconda::star \
    bioconda::kallisto \
    bioconda::bustools

mamba activate slsn
```

## Post-Installation Setup

### 1. Configuration
Edit the main configuration file:
```bash
cp config.yaml config_custom.yaml
# Edit config_custom.yaml with your specific paths and settings
```

### 2. Reference Genomes
You'll need to prepare reference genomes before running analyses. See the [Reference Genomes](5_reference_genomes.md) documentation for detailed instructions.

### 3. Sample Sheet
Create your sample sheet following the format described in [Sample Sheets](1_sample_sheets.md).

## Testing Your Installation

Run the pipeline on provided test data:
```bash
# Local run with 4 cores
snakemake -k -p --use-conda --conda-frontend mamba -j 4

# Test specific platform
snakemake out/test_sample/illumina/STAR/visium/counts_filtered/ -k -p --use-conda -j 4
```

## Hardware Requirements

### Minimum Requirements
- **CPU**: 4 cores
- **RAM**: 8GB
- **Storage**: 50GB free space

### Recommended Requirements
- **CPU**: 16+ cores
- **RAM**: 32GB+ (64GB for large datasets)
- **Storage**: 500GB+ SSD storage for optimal performance

### High-Performance Computing (HPC)
For large-scale analyses, consider using HPC clusters. `slide_snake` supports:
- SLURM job scheduler
- PBS/Torque
- SGE
- Cloud computing platforms (AWS, Google Cloud, Azure)

## Common Installation Issues

### Issue: Conda/Mamba solver conflicts
**Solution**: Use mamba instead of conda for faster dependency resolution:
```bash
conda install mamba -n base -c conda-forge
```

### Issue: Insufficient disk space
**Solution**: Configure temporary directory in config.yaml:
```yaml
TMPDIR: /path/to/large/temp/directory
```

### Issue: Memory errors during alignment
**Solution**: Adjust memory limits in config.yaml:
```yaml
MEMLIMIT_GB: 32G  # Adjust based on available RAM
```

## Next Steps

1. **Read the Documentation**: Start with [Sample Sheets](1_sample_sheets.md) to understand input requirements
2. **Explore Recipes**: Check [Recipes](2_recipes.md) to understand analysis options
3. **Prepare References**: Follow [Reference Genomes](5_reference_genomes.md) to set up your reference files
4. **Run Analysis**: Choose between [Short Read Pipeline](3a_short_read_pipeline.md) or [ONT Pipeline](3b_ont_pipeline.md)

## Getting Help

- **GitHub Issues**: Report bugs or request features at https://github.com/mckellardw/slide_snake/issues
- **Documentation**: All documentation files in the `docs/` directory
- **Test Data**: Use provided test datasets to validate your installation
