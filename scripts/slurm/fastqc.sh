#!/bin/bash
#SBATCH --job-name=fastqc
#SBATCH --time=0-3:00
#SBATCH --output=scripts/slurm/slurmOutput/R-%x.%j.log

### Arguments

#$1 = file containing input file paths (one per line)
#$2 = output directory

# set array correct for task
FILE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$1")

# run tool
fastqc "$FILE" -t ${SLURM_CPUS_PER_TASK} -o "$2"