#!/bin/bash
#SBATCH --job-name={{resources.job-name}}
#SBATCH --partition={{resources.partition}}
#SBATCH --time={{resources.time}}
#SBATCH --nodes={{resources.nodes}}
#SBATCH --cpus-per-task={{resources.threads}}
#SBATCH --mem={{resources.mem}}
#SBATCH --output={{resources.output}}

{{exec_job}}
