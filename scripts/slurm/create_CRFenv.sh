#!/bin/bash
#SBATCH --job-name=create_env
#SBATCH --time=01:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=2
#SBATCH --wait
#SBATCH --output=scripts/slurm/slurmOutput/R-%x.%j.log

conda env create -f /home/narosc/data/CoralReef_Fish/CRF_pipeline.yml