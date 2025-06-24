# Snakemake profiles for slide_snake

This directory contains Snakemake configuration profiles for running the slide_snake pipeline on different computing environments. Profiles allow you to customize resource allocation, scheduler settings, and execution parameters without modifying the main Snakefile.

## Available Profiles

- **`default/`** - Local execution with standard resource allocation
- **`slurm/`** - SLURM cluster execution with pe2 partition
- **`slurm_ne1/`** - SLURM cluster execution with cpu partition (alternative configuration)

## Using Profiles

To run the pipeline with a specific profile:

```bash
# Use default profile (local execution)
snakemake --profile profiles/default

# Use SLURM profile
snakemake --profile profiles/slurm

# Use alternative SLURM profile
snakemake --profile profiles/slurm_ne1
```

## Creating Custom Profiles

### 1. Profile Directory Structure

Create a new directory for your profile:
```
profiles/
├── your_profile_name/
│   ├── config.yaml      # Main configuration file
│   └── cluster-config.yaml  # Optional: cluster-specific settings
```

### 2. Configuration File Structure

The `config.yaml` file should contain the following sections:

#### A. Default Resources
Set baseline resource allocation for all rules:
```yaml
default-resources:
  runtime: "1:00:00"    # Maximum runtime (HH:MM:SS)
  mem_mb: 8192          # Memory in MB
  mem: "8G"             # Memory as string (for some schedulers)
  threads: 1            # Default thread count
  # Add scheduler-specific parameters as needed
  slurm_partition: "cpu"
  slurm_account: "your_account"
```

#### B. Thread Allocation
Customize thread usage for specific rules:
```yaml
set-threads:
  rule_name: 24          # Number of threads for this rule
  another_rule: 8
```

#### C. Resource Allocation
Override memory and other resources for specific rules:
```yaml
set-resources:
  rule_name:
    mem: "64G"
    mem_mb: 64000
    runtime: "4:00:00"
```

### 3. Rule Categories and Resource Guidelines

The slide_snake pipeline contains several categories of rules with different resource requirements:

#### Illumina Short-Read Processing (`ilmn_*`)
- **Preprocessing** (`ilmn_1*`): Light to moderate resources
  - FastQ merging: 16GB RAM, 8-24 threads
  - Adapter trimming: 16GB RAM, 8-24 threads
  - Quality control: 16GB RAM, 8-24 threads

- **rRNA Filtering** (`ilmn_2*`): Moderate resources
  - BWA alignment: 96GB RAM, 16 threads
  - Ribodetector: 32GB RAM, 16 threads

- **Genome Alignment** (`ilmn_3*`): High resources
  - STAR first pass: 128GB RAM, 24 threads
  - STAR second pass: 220-320GB RAM, 24 threads
  - Deduplication: 64GB RAM, 24 threads

- **Quantification** (`ilmn_4*`): High resources
  - Kallisto-bustools: 180-220GB RAM, 16 threads

- **Analysis** (`ilmn_5*-7*`): Variable resources
  - miRNA analysis: 200GB RAM, 4 threads
  - Data caching: 64GB RAM, 1 thread
  - Quality control: 8-16GB RAM, 8-56 threads

#### ONT Long-Read Processing (`ont_*`)
- **Preprocessing** (`ont_1a-1c`): Moderate to high resources
  - Adapter scanning: 64GB RAM, 56 threads
  - Read processing: 16-32GB RAM, 24 threads

- **Alignment** (`ont_1d-1f`): High resources
  - Minimap2 genome: 128GB RAM, 16 threads
  - ULTRA pipeline: 128GB RAM, 16 threads
  - Transcript quantification: 16-32GB RAM, 1-16 threads

- **Quality Control** (`ont_2*`): Light to moderate resources
  - Read QC: 8GB RAM, 16-24 threads
  - BAM QC: 32GB RAM, 1 thread

### 4. Environment-Specific Considerations

#### Local Execution
- Use conservative thread counts (≤ CPU cores)
- Monitor memory usage carefully
- Consider using smaller test datasets

#### SLURM Clusters
- Add SLURM-specific parameters:
```yaml
default-resources:
  slurm_partition: "cpu"        # or "gpu", "highmem", etc.
  slurm_account: "your_account"
  nodes: 1
  slurm_extra: "--constraint=avx2"  # Additional SLURM options
```

#### Other Schedulers (PBS, SGE, etc.)
- Replace `slurm_*` parameters with scheduler-specific ones
- Consult Snakemake documentation for your scheduler

### 5. Memory Considerations

The pipeline includes both `mem` (string) and `mem_mb` (integer) specifications for compatibility:
- Use `mem_mb` for precise memory allocation
- Use `mem` for scheduler-specific formatting
- Ensure both values are consistent

### 6. Runtime Estimation

Typical runtimes by rule category:
- Preprocessing: 30 minutes - 2 hours
- Alignment: 2-8 hours (depends on data size)
- Quantification: 1-4 hours
- Quality control: 30 minutes - 1 hour

Adjust `runtime` parameters based on your data size and system performance.

### 7. Testing Your Profile

1. **Start small**: Test with a subset of samples
2. **Monitor resources**: Check actual usage vs. allocated
3. **Adjust iteratively**: Fine-tune based on performance
4. **Document changes**: Keep notes on what works for your system

### 8. Example Custom Profile

Here's an example for a high-memory cluster:

```yaml
# profiles/highmem/config.yaml
default-resources:
  runtime: "2:00:00"
  mem_mb: 16384
  mem: "16G"
  threads: 1
  slurm_partition: "highmem"
  slurm_account: "my_project"

set-threads:
  ilmn_3a_STARsolo_secondPass: 32
  ont_1a_call_adapter_scan: 64
  
set-resources:
  ilmn_3a_STARsolo_secondPass:
    mem: "512G"
    mem_mb: 512000
    runtime: "12:00:00"
    
  ilmn_4a_kbpython_std:
    mem: "256G"  
    mem_mb: 256000
    runtime: "6:00:00"
```

### 9. Troubleshooting

Common issues and solutions:

- **Out of memory errors**: Increase `mem_mb` for failing rules
- **Job timeouts**: Increase `runtime` for long-running rules
- **Scheduler rejection**: Check partition limits and account permissions
- **Thread contention**: Reduce thread counts if system becomes unresponsive


## Contributing

If you create a useful profile for a common computing environment, consider contributing it back to the project via pull request.
