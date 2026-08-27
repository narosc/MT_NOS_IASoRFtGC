#!/bin/bash
#SBATCH --job-name=trim_galore
#SBATCH --time=0-3:00
#SBATCH --output=scripts/slurm/slurmOutput/R-%x.%j.log

### Arguments

#$1 = input file
#$2 = output folder
#$3 = read clipping 5'
#$4 = read clipping 3'
#$5 = adaptor contamination (1 (might overtrim), default 3 (moderate), 5+ (only if adaptor match strong))
    
# set array correct for task
FILE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$1")

# run tool
trim_galore "$FILE" -o "$2" --cores ${SLURM_CPUS_PER_TASK} \
    --clip_R1 "$3" --three_prime_clip_R1 "$4" --stringency "$5"