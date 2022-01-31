#!/bin/bash
#SBATCH --account=def-kagarwal
#SBATCH --mem-per-cpu=40G
#SBATCH --time=0-1:00
#SBATCH --array=0-124%25
#SBATCH --mail-user=simon.bernier@mail.mcgill.ca
#SBATCH --mail-type=ALL
runNumber=$(($SLURM_ARRAY_TASK_ID+125*5))
./gap-2dtfi $runNumber