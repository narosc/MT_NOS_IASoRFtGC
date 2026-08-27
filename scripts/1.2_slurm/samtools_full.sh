#!/bin/bash
#SBATCH --job-name=samtools_full
#SBATCH --time=02:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=4
#SBATCH --wait
#SBATCH --output=scripts/slurm/slurmOutput/R-%x.%j.log

### Arguments
# $1 = file containing input SAM paths
# $2 = output directory for bams and bai (index)
# $3 = output directory for flagstats and samtool stats

# set array correct for task
FILE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$1")

## extract sample name (keep environment clean and assure correct sample names)
sample=$(basename "$FILE" .sam)

## convert to BAM, sort (by genomic coordinates), mark duplicates and index
samtools view -@ "$SLURM_CPUS_PER_TASK" -b "$FILE" | \
  samtools sort -@ "$SLURM_CPUS_PER_TASK" -o "$2/${sample}.sorted.bam"

samtools markdup -@ "$SLURM_CPUS_PER_TASK" \
  "$2/${sample}.sorted.bam" \
  "$2/${sample}.sorted.markdup.bam"

samtools index "$2/${sample}.sorted.markdup.bam"

## run flagstat and SAMtools stats on BAM
samtools flagstat "$2/${sample}.sorted.markdup.bam" > "$3/${sample}.flagstat.txt"
samtools stats "$2/${sample}.sorted.markdup.bam" > "$3/${sample}.stats.txt"
samtools coverage "$2/${sample}.sorted.markdup.bam" > "$3/${sample}.coverage.txt"
samtools depth -aa "$2/${sample}.sorted.markdup.bam" > "$3/${sample}.depth.txt"
