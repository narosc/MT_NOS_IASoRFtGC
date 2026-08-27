## define path to -> query fasta file -> Zebrafish proteins
query_proteins <- 'intermed_outputs/danio_rerio_ref/Danio_rerio.GRCz11.pep.all.fa'

## define path to subject genome
subject_genome <- 'intermed_outputs/dascyllus_trimaculatus/ref/GCA_024666655.1.fna'

## create output directory
dir.create('intermed_outputs/dascyllus_trimaculatus/protein_annotation', recursive = TRUE, showWarnings = FALSE)

## build sbatch command
cmd <- paste(
  'sbatch ./scripts/slurm/blastn_prot_to_seq.sh',
  query_proteins,
  subject_genome,
  'intermed_outputs/dascyllus_trimaculatus/protein_annotation/drerio_vs_estriatus.tblastn.out'
)
system(cmd)
