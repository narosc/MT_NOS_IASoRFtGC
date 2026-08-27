#=============================================================================
### eighth (pre)analysis - Linkage decay estimation
#=============================================================================

# This script was written by Oliver M. Selmoni to estimate the linkage decay.
# This information is relevant before defining a bp window for mapping the SNPs
# to the genes.
# Inputs:
# 1) VCF file from fish species
# Output:
# 1. estimation of linkage decay

# Author: Oliver M. Selmoni (adjusted by Naroa Schweizer)
# last update: 20.08.26

# set correct working directory: setwd("C:/Users/Admin/Documents/Uni_Zuerich/Masters_Thesis/CoralReef_fish_local")

# load libraries
library(dartR)
library(vcfR)

#-----------------------------------------------------------------------------
## define paths and load files
#-----------------------------------------------------------------------------

# read vcf file
vcf  = read.vcfR('intermed_outputs/amphiprion_bicinctus/var/cohort_biallelic_fmiss_maf_pruned.vcf.gz')
# 'intermed_outputs/epinephelus_striatus/var/cohort_biallelic_fmiss_maf_pruned.vcf.gz'
# 'intermed_outputs/amphiprion_bicinctus/var/cohort_biallelic_fmiss_maf_pruned.vcf.gz'
# 'intermed_outputs/pterois_volitans/var/cohort_biallelic_fmiss_maf_pruned.vcf.gz'
# 'intermed_outputs/siphamia_tubifer/var/cohort_biallelic_fmiss_maf_pruned.vcf.gz'

# Convert to genlight format
gl = vcfR2genlight(vcf)

#=============================================================================
## 1) restructure genlight object
#=============================================================================

# recode SNPs number for the analysis (genelight doesn't like points in the identifiers)
gl$loc.names = paste0('snp',1:length(gl$loc.names))

# check number of SNPs per contig
table(gl$chromosome)

# retain only snps on the largest chromosome/contig
gl_filtered = gl.keep.loc(gl, loc.list = gl$loc.names[gl$chromosome=='NC_072778.1'])
# for Pterois volitans: 'CAMQGK010000052.1'
# for Amphiprion bicinctus: 'NC_072778.1'
# for Siphamia tubifer: 'CM035794.1'
# for Epinephelus striatus: 'CM069296.1'

# calculate LD between all SNPs on contig of interest (ncores = number of cpus available to parallelize)
LDreport = gl.report.ld(gl_filtered, ncores = 4)

# calculate distance (in nucleotides) between SNPs on contig of interest
LDreport$dist = abs(gl_filtered$position[LDreport$loc1]-gl_filtered$position[LDreport$loc2])

#=============================================================================
## 2) visualise analysis
#=============================================================================

# visualize how linkage (R2) changees over nucleotide distance
plot(LDreport$dist,
     LDreport$R2,
     pch=16,
     xlim=c(0,100000),
     cex=0.75,
     col=adjustcolor(1,0.3),
     ylab='Linkage (R2)',
     xlab='genomic distance (nucleotides)',
     main ='Linkage Decay for Amphiprion bicinctus')
### after ~10-20k nucleotides, SNPs seem to be mostly unlinked <- comment Oliver

#=============================================================================
## 3) quantify decay instead of eyeballing it (fit only, nothing added to the plot)
#=============================================================================

# fit an exponential decay curve to the raw pairs: R2 ~ a*exp(-b*dist) + c
# a = intercept at dist=0, c = background/noise floor, b = decay rate
max_dist  = 100000                                  # match the scatter plot's xlim
decay_fit = tryCatch(
  nls(R2 ~ a * exp(-b * dist) + c,
      data = LDreport[LDreport$dist > 0 & LDreport$dist <= max_dist, ],
      start = list(a = 0.5, b = 1/1000, c = 0.02),
      control = nls.control(maxiter = 200, warnOnly = TRUE)),
  error = function(e) NULL
)
# calculate and print LD decay rates
if (!is.null(decay_fit)) {
  
  co = coef(decay_fit)
  # distance at which the remaining excess LD is:
  # 50%, 25%, 10%, and 5% of its initial value
  decay_distances = c(
    '50% decay' = -log(0.50) / co['b'],
    '75% decay' = -log(0.25) / co['b'],
    '90% decay' = -log(0.10) / co['b'],
    '95% decay' = -log(0.05) / co['b']
  )
  # print only -> nothing from the fit is drawn on the plot above
  message('Distance (nucleotides) at which linkage decays by:')
  print(decay_distances)
  
}
