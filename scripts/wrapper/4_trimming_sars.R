## store read archive names in a vector
files <- list.files('intermed_outputs/ephinephelus_striatus/sra/', 
    pattern = "\\.fastq\\.gz$",
    full.names = TRUE
)

## create output directory
dir.create('intermed_outputs/ephinephelus_striatus/sra_trimmed/', showWarnings = FALSE) # showWarnings = FALSE -> no issue if dircetory already exists

## path to the file containing the list (same file used by the bash array)
list_file <- 'intermed_outputs/ephinephelus_striatus/sra_ids.txt'
writeLines(files, list_file)

## set trimming parameters:
c5 <- 10 # how many bp trimmed from 5'
c3 <- 1 # how many bp trimmed from 3'
s1 <- 1 # adaptor trimming stringency -> (1 default, 3 moderate, 5 stronger (5 bp overlap)) 

## define length of array
n <- length(files)

## command to trim all fastqc files according to clipping parameters
cmd <- paste(
    'sbatch',
    paste0('--array=1-', n, '%15'),
    '--cpus-per-task=4',
    '--mem=8G',
    './scripts/slurm/trim_galore.sh',
    list_file,
    'intermed_outputs/ephinephelus_striatus/sra_trimmed/',
    c5,
    c3,
    s1
)

system(cmd)
