#!/bin/bash
#SBATCH --job-name=bcf_mpileup
#SBATCH --time=0-3:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=16
#SBATCH --wait
#SBATCH --output=scripts/slurm/slurmOutput/R-%x.%j.log

### Arguments

# $1 = BAM list file (one BAM per line)
# $2 = output directory
# $3 = reference genome (.fna)
# $4 = max fraction missing per site (e.g. 0.1)
# $5 = MAF threshold (e.g. 0.05)

## joint mpileup, call variants and filter missing gt, biallelic and MAF (-Ou uncompressed vcf, -Oz compressed vcf gz)
# -f reference genome
# -b reference bam list for each sample
# -m2 -M2 together filter for only biallelic snps
# filter for MAF and missing gt -> -Oz -> vcf.gf compressed VCF 
# index file
# create statistics for quality of variant/SNP calling
bcftools mpileup \
  -f "$3" \
  -b "$1" \
  -Ou \
  --threads "$SLURM_CPUS_PER_TASK" |
bcftools call \
  -mv \
  -Ou |
bcftools view \
  -m2 -M2 -v snps \
  -Ou |
bcftools view \
  -i "F_MISSING<${4} && MAF>${5}" \
  -Oz \
  -o "$2/cohort_biallelic_fmiss_maf.vcf.gz"

bcftools index -t "$2/cohort_biallelic_fmiss_maf.vcf.gz"

bcftools stats \
  "$2/cohort_biallelic_fmiss_maf.vcf.gz" \
  > "$2/cohort_biallelic_fmiss_maf.stats.txt"
