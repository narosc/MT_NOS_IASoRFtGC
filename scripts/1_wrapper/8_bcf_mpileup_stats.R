## list sorted BAMs -> input
bams <- list.files('intermed_outputs/ephinephelus_striatus/aligned_bams/', full.names = TRUE, pattern='\\.sorted\\.bam$')  # exclude .sorted.bam.bai

## path for bam list file
bamlist_file <- 'intermed_outputs/ephinephelus_striatus/aligned_bams/bamlist.txt'

## write BAMs to text file, one per line
writeLines(bams, bamlist_file)

## create output directory
dir.create('intermed_outputs/ephinephelus_striatus/var', showWarnings = FALSE)

## define path to reference genome (.fna file)
refgen <- 'intermed_outputs/ephinephelus_striatus/ref/GCA_035609425.1.fna'

## define missing gt threshold
gt_miss_thresh <- 0.1

## define MAF threshold
MAF_thresh <- 0.05

## submit sbatch job
cmd <- paste(
    'sbatch ./scripts/slurm/bcf_mpileup.sh',
    shQuote(bamlist_file),
    'intermed_outputs/ephinephelus_striatus/var',
    refgen,
    gt_miss_thresh,
    MAF_thresh
)
system(cmd)