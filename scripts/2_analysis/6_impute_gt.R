#=============================================================================
### sixth analysis - impute genotype matrix for LFMM
#=============================================================================

# This script imputes the missing genotypes in the genlight matrix using the
# LEA package snmf function.
# Input:
# 1) genlight object filtered for missingnes, MAF and biallelic
# Output: 
# 1. an imputed genlight object with no missing values

# Notes: LEA::snmf() function is used to compute the statistically most likely 
# variant for NA values in samples (needs to know # K -> # populations if 
# structure is strong, missingness high). For my species so far this was not 
# the case.

# Author: Naroa Olivia Schweizer
# last update: 20.01.26

# set correct working directory: setwd("C:/Users/Admin/Documents/Uni_Zuerich/Masters_Thesis/CoralReef_fish_local")

# load libraries
library(LEA)            # adjust genotype matrix NULL values
library(adegenet)       # PCA & multivariate analysis of genetic data
library(ggplot2)        # general-purpose plotting

#-----------------------------------------------------------------------------
## define paths and load files
#-----------------------------------------------------------------------------

## filtered genligth object (biallelic and sites with less than 10% missing genotypes)
# define path to filtered and imputed genlight object
genlight_obj_path <- "intermed_outputs/epinephelus_striatus/var/estriatus_genlight_filtered.rds"

# read genlight object
genlight_obj <- readRDS(genlight_obj_path)

## outputs
# define output path
output_path <- "intermed_outputs/epinephelus_striatus/var"

#=============================================================================
## 1) impute missing values in genotype matrix
#=============================================================================

#-----------------------------------------------------------------------------
## 1.1) filter genlight object for NA values with LEA::snmf function
#-----------------------------------------------------------------------------

# load genlight object as df
gt <- as.data.frame(as.matrix(genlight_obj))

# make sure all values in matrix are numeric -> add row names (sample names) again
gt_num <- apply(gt, 2, function(col) as.numeric(col))
rownames(gt_num) <- rownames(gt)

# make NA values 9, for LEA::snmf() function
gt_num[is.na(gt_num)] <- 9

# lfmm formate of gt matrix (needed for LEA::snmf() function)
write.lfmm(gt_num, file.path(output_path, "gt.lfmm"))

# file path for snmf function
geno_file <- file.path(output_path, "gt.lfmm")

# compute snmf
project.snmf <- snmf(
  geno_file,
  K = 1:8,    # run for multiple K 
  ploidy = 2,
  entropy = TRUE,
  repetitions = 10,
  alpha = 100,
  project = "new"
)

#-----------------------------------------------------------------------------
## 1.2) visualise ancestral populations for species to select # K
#-----------------------------------------------------------------------------

# define range of K you actually ran
K_vals <- 1:max(project.snmf@K)

# extract minimal cross-entropy per K
ce_df <- data.frame(
  K = K_vals,
  cross_entropy = sapply(K_vals, function(k) {
    min(cross.entropy(project.snmf, K = k))
  })
)
# create plot for cross-entropy
ggplot(ce_df, aes(x = K, y = cross_entropy)) +
  geom_point(color = "grey40", size = 3) +
  geom_line(color = "grey40", linewidth = 1) +
  scale_x_continuous(breaks = ce_df$K) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Ancestral Population Cross-Entropy for Epinephelus striatus",
    x = "# ancestral populations [K]",
    y = "min cross-entropy"
  )
# NOTE: check where cross-entropy is minimal -> then use this for # K

#-----------------------------------------------------------------------------
## 1.3) impute using K with min cross entropy
#-----------------------------------------------------------------------------

# select the run with the lowest cross-entropy value
K_best <- ce_df$K[which.min(ce_df$cross_entropy)]
# Pterois volitans: K = 2 
# Amphiprion bicinctus: K = 3
# Siphamia tubifer: K = 1
# Epinephelus striatus: K = 1
best_run <- which.min(cross.entropy(project.snmf, K = K_best))

# impute the missing genotypes
impute(
  project.snmf,
  geno_file,
  method = "mode",
  K = K_best,
  run = best_run
)

# define path to imputed file
imputed_file <- file.path(output_path, "gt.lfmm_imputed.lfmm")

# read imputed genotype matrix
imputed_gt <- as.matrix(read.table(imputed_file, header = FALSE))

# add columns and rows back from original gt (identify sample and contig)
rownames(imputed_gt) <- rownames(gt_num)
colnames(imputed_gt) <- colnames(gt_num)

# rebuild properly
genlight_fixed <- as.genlight(imputed_gt)

# set ploidy (diploid species)
ploidy(genlight_fixed) <- 2

# restore metadata
indNames(genlight_fixed) <- rownames(imputed_gt)
locNames(genlight_fixed) <- colnames(imputed_gt)

# restore population information
pop(genlight_fixed) <- pop(genlight_obj)

# if chromosome/position existed, restore them too
chromosome(genlight_fixed) <- chromosome(genlight_obj)
position(genlight_fixed) <- position(genlight_obj)

#=============================================================================
## 2) save outputs
#=============================================================================

# save imputed and filtered genlight object
saveRDS(genlight_fixed,
        file = file.path(output_path, "epistri_genlight_filtered_imputed.rds"))
