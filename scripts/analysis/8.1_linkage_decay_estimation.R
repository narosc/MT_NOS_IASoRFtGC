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
# last update: 12.05.26

# set correct working directory: setwd("C:/Users/Admin/Documents/Uni_Zuerich/Masters_Thesis/CoralReef_fish_local")

# load libraries
library(dartR)
library(vcfR)

#-----------------------------------------------------------------------------
## define paths and load files
#-----------------------------------------------------------------------------

### Read vcf file
vcf  = read.vcfR('intermed_outputs/epinephelus_striatus/var/cohort_biallelic_fmiss_maf.vcf.gz')

### Convert to genlight format
gl = vcfR2genlight(vcf)

#=============================================================================
## 1) restructure genlight object
#=============================================================================

### Recode SNPs number for the analysis (genelight doesn't like points in the identifiers)
gl$loc.names = paste0('snp',1:length(gl$loc.names))

## check number of SNPs per contig
table(gl$chromosome)

### retain only snps on the largest chromosome/contig
gl_filtered = gl.keep.loc(gl, loc.list = gl$loc.names[gl$chromosome=='CM069296.1'])
# for Pterois volitans: 'CAMQGK010000052.1'
# for Amphiprion bicinctus: 'NC_072778.1' -> 'NC_072766.1'
# for Siphamia tubifer: 'CM035794.1'
# for Epinephelus striatus: 'CM069296.1'

### calculate LD between all SNPs on contig of interest (ncores = number of cpus available to parallelize)
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
     main ='Linkage Decay for Epinephelus striatus')
### after ~10-20k nucleotides, SNPs seem to be mostly unlinked 

#=============================================================================
## 3) quantify decay instead of eyeballing it (bin + fit on the same LDreport)
#=============================================================================

# bin R2 by distance and take the mean per bin - the raw scatter above is
# heavily overplotted near dist=0, so a binned mean makes the decay trend
# visible instead of having to read it off point density by eye
bin_width = 500                                     # nucleotides per bin
max_dist  = 100000                                  # match the scatter plot's xlim

LDreport_bin = LDreport[LDreport$dist <= max_dist, ]
LDreport_bin$bin = cut(LDreport_bin$dist, breaks = seq(0, max_dist, by = bin_width), include.lowest = TRUE)

bin_means = aggregate(R2 ~ bin, data = LDreport_bin, FUN = mean)
bin_means$dist_mid = seq(bin_width / 2, max_dist - bin_width / 2, by = bin_width)[seq_len(nrow(bin_means))]

# same raw scatter as section 2, but lighter, with the binned mean R2 overlaid
plot(LDreport$dist,
     LDreport$R2,
     pch=16,
     xlim=c(0,100000),
     cex=0.75,
     col=adjustcolor(1,0.15),
     ylab='Linkage (R2)',
     xlab='genomic distance (nucleotides)',
     main ='Linkage Decay for Epinephelus striatus (binned mean + fitted decay)')
lines(bin_means$dist_mid, bin_means$R2, col = 'red', lwd = 2)

# fit an exponential decay curve to the raw (unbinned) pairs:
# R2 ~ a*exp(-b*dist) + c
# a+c = intercept at dist=0, c = background/noise floor, b = decay rate
decay_fit = tryCatch(
  nls(R2 ~ a * exp(-b * dist) + c,
      data = LDreport[LDreport$dist > 0 & LDreport$dist <= max_dist, ],
      start = list(a = 0.5, b = 1/2000, c = 0.02),
      control = nls.control(maxiter = 200, warnOnly = TRUE)),
  error = function(e) NULL
)

if (!is.null(decay_fit)) {
  
  co = coef(decay_fit)
  half_decay_dist = -log(0.5) / co['b']    # distance at which excess LD (above background) halves
  
  fit_x = seq(1, max_dist, by = 100)
  lines(fit_x, predict(decay_fit, newdata = data.frame(dist = fit_x)), col = 'blue', lwd = 2, lty = 2)
  
  legend('topright',
         legend = c('binned mean R2', paste0('fitted decay curve (half-decay = ', round(half_decay_dist), ' bp)')),
         col = c('red', 'blue'), lwd = 2, lty = c(1, 2), bty = 'n')
  
  half_decay_dist   # <- report this number instead of an eyeballed "±X kb"
  
} else {
  warning('exponential decay fit did not converge - report the binned mean curve visually instead of a half-decay distance for this species')
}

### NOTE: this fit assumes R2 approaches a stable background level c as distance
### increases. For a species/chromosome with elevated, non-decaying LD across
### the whole window (see LD write-up for Amphiprion bicinctus), the model may
### not converge, or may return an unstable/very large half_decay_dist - if so,
### that itself is worth reporting (i.e. "no reliable exponential decay could be
### fitted"), rather than forcing a number that doesn't reflect the data.

