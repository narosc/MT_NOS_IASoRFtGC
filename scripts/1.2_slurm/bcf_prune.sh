#!/bin/bash
#SBATCH --job-name=bcf_prune
#SBATCH --time=01:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4
#SBATCH --wait
#SBATCH --output=scripts/slurm/slurmOutput/R-%x.%j.log

# $1 = input vcf.gz
# $2 = output vcf.gz
# $3 = r2 threshold -> from 0 to 1
# $4 = window size (in bp)

bcftools +prune \
  "$1" \
  -m "$3" \
  -w "$4" \
  -Oz \
  -o "$2"

bcftools index -t "$2"