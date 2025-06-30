# Installing `slide_snake`

`slide_snake` is a flexible Snakemake workflow for processing spatial RNA-sequencing data from various commercial and custom platforms. This guide will walk you through the installation process and initial setup.

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
snakemake --use-conda --conda-frontend mamba -j 4
```
- *Note* - the first run will take much longer as all of the rule-specific conda environments need to be built.

## Post-Installation Setup

### 1. Configuration
Edit the main configuration file:
- Update `SAMPLE_SHEET_PATH` to point to your actual sample sheet (default will point to the included toy datasets)

### 2. Reference Genomes
You'll need to prepare `STAR` and/or `kallisto` references before running analyses. See the [Reference Genomes](5_reference_genomes.md) documentation for detailed instructions.

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

### Recommended Requirements
- **CPU**: 16+ cores
- **RAM**: 32GB+ (64GB for large datasets)

## Common Installation Issues

### Issue: Memory errors during alignment
**Solution**: Adjust memory limits in config.yaml:
```yaml
MEMLIMIT_GB: 32G  # Adjust based on available RAM
```

## Next Steps

1. **Read the Documentation**: Start with [Sample Sheets](1_sample_sheets.md) to understand input requirements
2. **Explore Recipes**: Check [Recipes](2a_recipes.md) to understand analysis options
3. **Prepare References**: Follow [Reference Genomes](5_reference_genomes.md) to set up your reference files
4. **Run Analysis**: Choose between [Short Read Pipeline](3a_short_read_pipeline.md) or [ONT Pipeline](3b_ont_pipeline.md)

## Getting Help

- **GitHub Issues**: Report bugs or request features at https://github.com/mckellardw/slide_snake/issues
- **Documentation**: All documentation files in the `docs/` directory
- **Test Data**: Use provided test datasets to validate your installation
