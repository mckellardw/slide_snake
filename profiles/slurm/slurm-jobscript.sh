#!/bin/bash
#SBATCH --account={{resources.account}}
#SBATCH --time={{resources.time}}
#SBATCH --partition={{resources.partition}}
{{exec_job}}
