#!/bin/bash
#SBATCH --job-name=parallel_fastq_dump
#SBATCH --time=0-8:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=10
#SBATCH --array=1-281%20
#SBATCH --output=scripts/slurm/slurmOutput/R-%x.%j.log

# notes: array 281 (samples in total for 3 species) % 20 (max 20 jobs in parallel)

### Arguments

#$1 = sra ID
#$2 = output directory

# get the correct SRA ID for this array task
SRA_ID=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$1")

# skip sample if already downloaded
if [ -f "$2/${SRA_ID}.fastq.gz" ]; then
    exit 0
fi

# run parallel-fastq-dump
parallel-fastq-dump \
  --sra-id "${SRA_ID}" \
  --threads ${SLURM_CPUS_PER_TASK} \
  --outdir "$2" \
  --gzip \
  --clip \
  --skip-technical