#!/bin/bash
#SBATCH --job-name=read_align_ngm
#SBATCH --time=01:00:00
#SBATCH --output=scripts/slurm/slurmOutput/R-%x.%j.log

### Arguments

# $1 = path to reference genome
# $2 = file containing input file paths
# $3 = output directory

# set array correct for task
FILE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$2")

# derive sample name
base=$(basename "$FILE")
base=${base%%_trimmed*}              # remove "_trimmed"
sample_name=${base%.fq.gz}           # remove .fq.gz

# define output SAM path
OUT="$3/${sample_name}.sam"

# run tool
ngm -r "$1" -q "$FILE" -o "$OUT" -t ${SLURM_CPUS_PER_TASK}