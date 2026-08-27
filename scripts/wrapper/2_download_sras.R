## load metadata of bioproject
sra_table <- read.csv('data/SraRunTable_PRJEB36904_Ephinephelusstriatus.csv', 
                    stringsAsFactors = FALSE)

## create output directory
dir.create('intermed_outputs/ephinephelus_striatus/sra', showWarnings = FALSE) # showWarnings = FALSE -> no issue if dircetory already exists

## extract sra IDs -> SRRXXXXX
sraIDs <- sra_table$Run

## set a subset to test if working
# sraIDs <- sraIDs[c(1, 30, 60, 100)] # set different ones, to see if overall samples similar quality

## write SRA IDs to file (needed for SLURM array)
writeLines(
  sraIDs,
  'intermed_outputs/ephinephelus_striatus/sra_ids.txt'
)

## count number of samples (for array range)
n_samples <- length(sraIDs)

## command to download sra samples
cmd <- paste(
  'sbatch',
  paste0('--array=1-', n_samples, '%20'),
  './scripts/slurm/parallel_fastq_dump.sh',
  'intermed_outputs/ephinephelus_striatus/sra_ids.txt',
  'intermed_outputs/ephinephelus_striatus/sra/'
)

system(cmd)