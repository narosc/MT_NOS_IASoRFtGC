## create output directory
dir.create('intermed_outputs/dascyllus_trimaculatus/ref', showWarnings = FALSE) # showWarnings = FALSE -> no issue if dircetory already exists

## define which genbank assembly accession number
gca <- 'GCA_024666655.1' # Dascyllus trimaculatus as reference for species 5 (Dascyllus trimaculatus)

## download reference genome from ncbi database
cmd <- paste(
  'sbatch ./scripts/slurm/download_ref_genome.sh',
  gca,
  'intermed_outputs/dascyllus_trimaculatus/ref'
  )
  system(cmd) # submitt job

  ## list of reference genome genbank assembly numbers
  # 'GCA_947000775.1' # pterois miles as reference for species 1 (Pterois volitans)
  # 'GCF_022539595.1' # amphiprion ocellaris as reference for species 2 (Amphiprion bicinctus)
  # 'GCA_020466265.1' # siphamia tubifer as reference for species 3 (Siphamia tubifer)
  # 'GCA_035609425.1' # Epinephelus awoara as reference for species 4 (Epinephelus striatus)
  # 'GCA_024666655.1' # Dascyllus trimaculatus as reference for species 5 (Dascyllus trimaculatus)