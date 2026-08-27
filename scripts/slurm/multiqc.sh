#!/bin/bash
#SBATCH --job-name=multiqc
#SBATCH --time=00:30:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=1
#SBATCH --wait
#SBATCH --output=scripts/slurm/slurmOutput/R-%x.%j.log

### Arguments

# $1 = input directory (containing fastqc, flagstats etc. -> intermed_outputs)
# $2 = output directory

# run MultiQC on the input directory
multiqc "$1" -o "$2"