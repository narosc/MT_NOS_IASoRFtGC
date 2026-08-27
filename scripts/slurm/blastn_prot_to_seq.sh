#!/bin/bash
#SBATCH --job-name=blast_prot_to_seq
#SBATCH --time=48:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=16
#SBATCH --wait
#SBATCH --output=scripts/slurm/slurmOutput/R-%x.%j.log

### Arguments
# $1 = protein fasta (query)
# $2 = subject genome (SOI)
# $3 = output file

# define DB prefix -> strip .fna suffix
DB_PREFIX="${2%.fna}"

# build database (only if it does not exist)
if [ ! -f "${DB_PREFIX}.nin" ]; then
    makeblastdb \
        -in "$2" \
        -dbtype nucl \
        -out "$DB_PREFIX"
fi

# run tblastn
tblastn \
  -query "$1" \
  -db "$DB_PREFIX" \
  -out "$3" \
  -evalue 5e-07 \
  -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore" \
  -num_threads "$SLURM_CPUS_PER_TASK"