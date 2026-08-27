#=============================================================================
### nineth analysis - PicMin - across-species analysis
#=============================================================================

# This script runs the across-species analysis using PicMin package. It looks
# for a signal repeated signal across species linking genetic variation to 
# environmental variable.
# Inputs: 
# 1) map_SNPs_genes outputs per species for the specific GEA method (LFMM or CA)
# 2) costom PicMin function written by Oliver Selmoni
# Output: 
# 1. picmin results RData list per environmental variable

# set correct working directory: setwd("C:/Users/Admin/Documents/Uni_Zuerich/Masters_Thesis/CoralReef_fish_local")

# Author: Naroa Olivia Schweizer
# last update: 02.07.26

# load libraries 
library(PicMin)          # runs PicMin function
library(foreach)         # needed to create forloop 
library(doParallel)      # needed to run picmin in parallel
library(dplyr)           # data wrangling (filter, summarize, join)
library(tidyverse)      # restructure tables and formats

#-----------------------------------------------------------------------------
## define paths and load files
#-----------------------------------------------------------------------------

# Note: one has to specify which version to use (_LFMM.RData or _CA.RData)

## R-Object of combined p-values per genes
# Pterois volitans
# define path to r-object
gene_path <- "intermed_outputs/pterois_volitans/protein_annotation/SNPs_mapped_to_genes_LFMM.RData"

# load r-object (gene_results)
load(gene_path)
ptevol_gene_results <- gene_results

## Amphiprion bicinctus
# define path to r-object
gene_path <- "intermed_outputs/amphiprion_bicinctus/protein_annotation/SNPs_mapped_to_genes_LFMM.RData"

# load r-object (gene_results)
load(gene_path)
amphibic_gene_results <- gene_results

## Siphamia tubifer
# define path to r-object
gene_path <- "intermed_outputs/siphamia_tubifer/protein_annotation/SNPs_mapped_to_genes_LFMM.RData"

# load r-object (gene_results)
load(gene_path)
siphtub_gene_results <- gene_results

## Epinephelus striatus
# define path to r-object
gene_path <- "intermed_outputs/epinephelus_striatus/protein_annotation/SNPs_mapped_to_genes_LFMM.RData"

# load r-object (gene_results)
load(gene_path)
epistri_gene_results <- gene_results

# create vectors containing prefix names of all species
prefix_fish <- c("ptevol","amphibic", "siphtub", "epistri")

## outputs
# create output direcotry for lfmm outputs
# dir.create('intermed_outputs/picmin/LFMM', showWarnings = FALSE)
# define output path to save etc.
output_path <- "intermed_outputs/picmin/LFMM"

# load custom function 
source('scripts/analysis/9.1_custom_picmin_functions_OS.R')

#=============================================================================
## 1) prepare df as input for PicMin
#=============================================================================

#-----------------------------------------------------------------------------
## 1.1) create df for each env variable with combined p per species
#-----------------------------------------------------------------------------

# combine all four species' Rdata lists into one, named by species prefix
picmin_inputs_all <- list(
  ptevol   = ptevol_gene_results,
  amphibic = amphibic_gene_results,
  siphtub  = siphtub_gene_results,
  epistri  = epistri_gene_results
)

# create function to extract df per env variable containing emp p-values of each species
merge_by_environment <- function(input_list) {
  # get shared environment names across all species
  env_names <- Reduce(intersect, lapply(input_list, names))
  species_names <- names(input_list)
  result <- lapply(env_names, function(env) {
    df_list <- lapply(species_names, function(species) {
      prefix <- species  # already named by prefix in picmin_inputs_all
      
      df <- input_list[[species]][[env]] %>%
        dplyr::filter(gene != "") |>              # drop unassigned SNPs
        dplyr::select(gene, combined_p) %>%
        dplyr::group_by(gene) %>%
        dplyr::summarise(combined_p = min(combined_p, na.rm = TRUE), .groups = "drop") %>%
        dplyr::distinct()
      colnames(df)[colnames(df) == "combined_p"] <- prefix
      return(df)
    })
    # merge across species by gene
    merged_df <- Reduce(function(x, y) {
      dplyr::full_join(x, y, by = "gene")
    }, df_list)
    return(merged_df)
  })
  names(result) <- env_names
  return(result)
}

# run function across all species/environments
merged_env_picmin_input <- merge_by_environment(picmin_inputs_all)

#-----------------------------------------------------------------------------
## 1.3) save output list
#-----------------------------------------------------------------------------

## save merged_env_dfs list as an r-object list
# define path and name!
output_file1 <- file.path(
  output_path,
  "ptevol_amphibic_siphtub_epistri_picmin_input_combinedP_LFMM.RData"
)
# save
save(merged_env_picmin_input, file = output_file1)

#=============================================================================
## 2) run PicMin function
#=============================================================================

#-----------------------------------------------------------------------------
## 2.1) combine all p-values
#-----------------------------------------------------------------------------

# ## if already calculated just load this and function
# # load picmin input data
# load('intermed_outputs/picmin/CA/ptevol_amphibic_siphtub_picmin_input_combinedP_CA.RData')
# load('intermed_outputs/picmin/LFMM/ptevol_amphibic_siphtub_picmin_input_combinedP_lfmm.RData')

## merge pvalues in a list
PVALS_ALL = merged_env_picmin_input
# make gene_id be rownames not a seperate column
PVALS_ALL <- lapply(merged_env_picmin_input, function(df) {
  df |>
    as.data.frame() |>
    tibble::column_to_rownames(var = "gene")
})

#-----------------------------------------------------------------------------
## 2.2) set up parallel backend for processor useage
#-----------------------------------------------------------------------------

#setup parallel backend to use many processors
# to see how many cores: parallel::detectCores()
# [1] 8 (naroa)
cores_to_use = 8
cl <- makeCluster(cores_to_use) 
registerDoParallel(cl)

#-----------------------------------------------------------------------------
## 2.3) run PicMin with function from Oliver
#-----------------------------------------------------------------------------

# record starting time
t0 = Sys.time()

# precompute null distribution of p-values for 4 species
nullP <- PicMinNull(linMin = 3, linMax = 4)

# create empty results list, later store picmin output per environment
picmin_results <- list()

## run PicMin for every environmental variable
for (env in names(PVALS_ALL)) {
  
  # p-values table for current environmental variable
  # (gene IDs are already rownames from the PVALS_ALL creation step above --
  #  no need to re-set rownames or drop a column here)
  all_lins_p <- PVALS_ALL[[env]]
  
  # convert every species column to empirical p-values
  for (col_i in seq_len(ncol(all_lins_p))) {
    all_lins_p[, col_i] <- PicMin:::EmpiricalPs(all_lins_p[, col_i], large_i_small_p = FALSE)
  }
  
  ### Run picmin
  my_picmin_output <- RunPicmin(all_lins_p, nullP = nullP)
  
  print(env)
  print(sum(my_picmin_output$pooled_q < 0.5, na.rm = TRUE))
  
  # store directly in the results list -- no per-environment save
  picmin_results[[env]] <- my_picmin_output
}

# stop cluster once all environments are done
stopCluster(cl)

# report runtime
print(Sys.time() - t0)
# for both CA and LFMM between 8-9 minutes

#-----------------------------------------------------------------------------
## 2.4) save picmin results
#-----------------------------------------------------------------------------

# save the whole results object once, at the end
save(picmin_results, file = "picmin_results_lfmm.rda")

#-----------------------------------------------------------------------------
## 2.5) check output
#-----------------------------------------------------------------------------

# check for genes with q-value below 0.5
for (name in names(picmin_results)) {
  df <- picmin_results[[name]]
  
  sig <- df[df$q < 0.5, c("locus", "p", "q", "n_est")]
  
  if (nrow(sig) > 0) {
    cat("\n###", name, "###\n")
    print(sig)
  }
}
