## define input directory
input_dir <- 'intermed_outputs/ephinephelus_striatus/' # will search through all subfolders

## create output directory
dir.create('intermed_outputs/ephinephelus_striatus/multiqc_raw', showWarnings = FALSE)

## build sbatch command
cmd <- paste(
  'sbatch ./scripts/slurm/multiqc.sh',
  input_dir,
  'intermed_outputs/ephinephelus_striatus/multiqc_raw'
)
system(cmd) # submitt job