# envs

Environments were built and exported using:
```
conda export env --no-builds > envs/env.yml
```
- The `-no-builds` flag is important to make the envs work cross-platform
- Include the `--use-conda` flag in your snakemake call to use these
- If you prefer mamba, you can use the `--