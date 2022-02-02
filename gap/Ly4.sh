#!/bin/bash
#SBATCH --account=def-kagarwal
#SBATCH --mem-per-cpu=8G
#SBATCH --time=0-1:00
#SBATCH --array=0-124%25
#SBATCH --mail-user=simon.bernier@mail.mcgill.ca
#SBATCH --mail-type=ALL
runNumber=$(($SLURM_ARRAY_TASK_ID))
./gap-2dtfi $runNumber