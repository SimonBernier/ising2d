#!/bin/bash
#SBATCH --account=def-kagarwal
#SBATCH --mem-per-cpu=1G
#SBATCH --time=0-1:00
#SBATCH --array=0-57
#SBATCH --mail-user=simon.bernier@mail.mcgill.ca
#SBATCH --mail-type=ALL
runNumber=$$(($SLURM_ARRAY_TASK_ID+58*2))
./ising2-PT $runNumber