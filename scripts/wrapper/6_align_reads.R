## define paths
ref_genome <- 'intermed_outputs/ephinephelus_striatus/ref/GCA_035609425.1.fna'
files <- list.files(
  'intermed_outputs/ephinephelus_striatus/sra_trimmed/', 
  pattern = "\\.fq\\.gz$", 
  full.names = TRUE
)

## create output directory
dir.create('intermed_outputs/ephinephelus_striatus/aligned_reads', showWarnings = FALSE)

## path to the file containing the list (same file used by the bash array)
list_file <- 'intermed_outputs/ephinephelus_striatus/sra_ids.txt'
writeLines(files, list_file)

## define length of array
n <- length(files)

## command to align each cleand read to the reference genome
cmd <- paste(
  'sbatch',
  paste0('--array=1-', n, '%15'), 
  '--cpus-per-task=8',
  '--mem=32G',
  './scripts/slurm/read_alignment_ngm.sh',
  shQuote(ref_genome),
  shQuote(list_file),
  shQuote('intermed_outputs/ephinephelus_striatus/aligned_reads')
)

system(cmd)
