#!/bin/bash
#SBATCH --job-name=slide_snake                 # Job name
#SBATCH --partition=pe2                     # Partition Name
#SBATCH --mail-type=ALL                 # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=ecordina@nygenome.org      # Where to send mail	
#SBATCH --mem=16G                            # Job memory request. Different units can be specified using the suffix [K|M|G|T]
#SBATCH --time=15:00:00                       # Time limit 4 hours
#SBATCH --output=stdout_%j.log               # Standard output and error log
#SBATCH --cpus-per-task=4


cd /gpfs/commons/home/ecordina/uST/github/slide_snake/
snakemake --cluster-config config/slurm.yaml --cluster "sbatch -p {cluster.partition} -t {cluster.time} -N {cluster.nodes} --mem {cluster.mem} -o {cluster.output} --cpus-per-task={cluster.threads}" -j 32 --cluster-cancel scancel --use-conda --conda-frontend mamba
