## pruning parameters (r2 between 0 to 1, "kernel"/window in bp)
r2_threshold <- 0.9
window_size <- 1000

## submit sbatch job
cmd <- paste(
  "sbatch",
  "./scripts/slurm/bcf_prune.sh",
  shQuote("intermed_outputs/ephinephelus_striatus/var/cohort_biallelic_fmiss_maf.vcf.gz"),
  shQuote("intermed_outputs/ephinephelus_striatus/var/cohort_biallelic_fmiss_maf_pruned.vcf.gz"),
  r2_threshold,
  window_size
)

system(cmd)