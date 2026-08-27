## NOTE: this might not be needed for all species... only if some reads didn't download

## load metadata of bioproject
sra_table <- read.csv('data/SraRunTable_PRJNA385011_Siphamiatubifer.csv', 
                    stringsAsFactors = FALSE)

## create output directory
dir.create('intermed_outputs/siphamia_tubifer/sra', showWarnings = FALSE) # showWarnings = FALSE -> no issue if dircetory already exists

## extract sra IDs -> SRRXXXXX
sraIDs <- sra_table$Run

## missing SRAs
sraIDs_missing <- c(
 "SRR5489652",
 "SRR5489724",
 "SRR5489773"
)

## write SRA IDs to file (needed for SLURM array)
writeLines(
  sraIDs_missing,
  "intermed_outputs/siphamia_tubifer/sra/sra_ids_miss.txt"
)

## count number of samples (for array range)
n_samples <- length(sraIDs_missing)

## command to download sra samples
cmd <- paste(
  "sbatch",
  paste0("--array=1-", n_samples),
  "./scripts/slurm/parallel_fastq_dump.sh",
  "intermed_outputs/siphamia_tubifer/sra/sra_ids_miss.txt",
  "intermed_outputs/siphamia_tubifer/sra/"
)

system(cmd)