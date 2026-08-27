## list SAMs
files <- list.files(
  'intermed_outputs/ephinephelus_striatus/aligned_reads/', 
  pattern = "\\.sam$", 
  full.names = TRUE
)

## create two output directories (one for the bams/bais and one for the stats)
dir.create('intermed_outputs/ephinephelus_striatus/aligned_bams', showWarnings = FALSE)
dir.create('intermed_outputs/ephinephelus_striatus/samtools_stats', showWarnings = FALSE)

## path to the file containing the list (same file used by the bash array)
list_file <- 'intermed_outputs/ephinephelus_striatus/sra_ids.txt'
writeLines(files, list_file)

## define length of array
n <- length(files)

## command to get samtools stats and flagstats for all bam files
cmd <- paste(
  'sbatch',
  paste0('--array=1-', n, '%15'),  # 15 tasks max at once
  '--cpus-per-task=4',
  '--mem=8G',
  './scripts/slurm/samtools_full.sh',
  shQuote(list_file),
  shQuote('intermed_outputs/ephinephelus_striatus/aligned_bams'),
  shQuote('intermed_outputs/ephinephelus_striatus/samtools_stats')
)

system(cmd)