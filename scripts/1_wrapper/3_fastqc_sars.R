## store read archive names in a vector
files <- list.files('intermed_outputs/ephinephelus_striatus/sra/', 
    pattern = "\\.fastq\\.gz$",
    full.names = TRUE
)

## create output directory
dir.create('intermed_outputs/ephinephelus_striatus/fastqc', showWarnings = FALSE) # showWarnings = FALSE -> no issue if dircetory already exists

## path to the file containing the list (same file used by the bash array)
list_file <- 'intermed_outputs/ephinephelus_striatus/sra_ids.txt'
writeLines(files, list_file)

## define length of array
n <- length(files)

## command to run fastqc for each file
cmd <- paste(
  'sbatch',
  paste0('--array=1-', n, '%15'),
  '--cpus-per-task=4',
  '--mem=4G',
  './scripts/slurm/fastqc.sh',
  list_file,
  'intermed_outputs/ephinephelus_striatus/fastqc'
)

system(cmd)