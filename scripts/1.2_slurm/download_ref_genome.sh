#!/bin/bash
#SBATCH --job-name=download_ref
#SBATCH --time=01:00:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=2
#SBATCH --wait
#SBATCH --output=scripts/slurm/slurmOutput/R-%x.%j.log

### Arguments

#$1 = input file
#$2 = output directory

## create output directory with name specified in wrapper script
mkdir -p "$2"

## download the reference genome zip to the correct location
datasets download genome accession "$1" --include genome --filename "$2/ncbi_dataset.zip"

## unzip to a subfolder to access fna file
unzip "$2/ncbi_dataset.zip" -d "$2/ncbi_dataset/"

## find and capture FASTA file inside unziped folder, then rename and move it to data/ref
fna_file=$(find "$2/ncbi_dataset" -type f -name "*.fna")
mv "$fna_file" "$2/$1.fna"

## remove the zip file -> not necessairy anymore
rm "$2/ncbi_dataset.zip"
rm -r "$2/ncbi_dataset"