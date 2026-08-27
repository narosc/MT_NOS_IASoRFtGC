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

##----
#old
#-----

# #-----------------------------------------------------------------------------
# ## 2.4) stopp cluster
# #-----------------------------------------------------------------------------
# 
# ### Stop Cluster  
# stopCluster(cl)
# 
# #=============================================================================
# ## 3) save PicMin results
# #=============================================================================
# 
# #-----------------------------------------------------------------------------
# ## 3.1) create output list by environment and save results
# #-----------------------------------------------------------------------------
# 
# # create results list by environment
# picmin_results <- list()
# 
# for (env in names(PVALS_ALL)) {
#   
#   load(paste0(env, "_my_picmin_output.rda"))  # loads object 'my_picmin_output'
#   
#   picmin_results[[env]] <- my_picmin_output
# }
# 
# # save output
# save(picmin_results, file = "picmin_results_ca.rda")
# 
# #=============================================================================
# ## 4) check outputs
# #=============================================================================
# 
# # if analysis already ran and results are saved -> reload them:
# # define path to r-object of picmin results
# picmin_path <- "picmin_results_lfmm.rda"
# 
# # load r-object (picmin_path)
# load(picmin_path)
# 
# #-----------------------------------------------------------------------------
# ## 4.1) check results with specific q-value
# #-----------------------------------------------------------------------------
# 
# # check for genes with q-value below 0.5
# for (name in names(picmin_results)) {
#   df <- picmin_results[[name]]
#   
#   sig <- df[df$q < 0.5, c("locus", "p", "q", "n_est")]
#   
#   if (nrow(sig) > 0) {
#     cat("\n###", name, "###\n")
#     print(sig)
#   }
# }
# 
# # # check output
# # for (env in names(PVALS_ALL)) { # for every environmental variable...
# #  
# #   load(paste0(env,'_my_picmin_output.rda'))
# #   print(env)
# #   print(my_picmin_output[my_picmin_output$pooled_q<0.3,])
# #   
# #    
# # }
# # 
# # pvals = merged_env_dfs[['s_allDB_URBA_avg']]
# # rownames(pvals) = pvals$gene
# # pvals['ENSDARG00000098183',]
# # 
# # 
# # ### Compare with my results
# # load('env1_picMin_results.rda')
# # load('env1_my_picmin_output.rda')
# # 
# # head(picMin_results)
# # head(my_picmin_output)
# # 
# # ## check if same order of p-values
# # rownames(picMin_results) = picMin_results$locus
# # rownames(my_picmin_output) = my_picmin_output$locus
# # 
# # ## check correlation p-values
# # plot(picMin_results$p, my_picmin_output[picMin_results$locus,'p'], pch=16, cex=0.3)
# # 
# # ## Check manhattan plots
# # plot(-log(picMin_results$p))
# # plot(-log(my_picmin_output$p))
# 
reef_vibrant <- c(
  "#154C79",  # Midnight Blue
  "#2D7DD2",  # Reef Blue
  "#1FC8C8",  # Tropical Cyan
  "#78E3D8",  # Seafoam
  "#F46D75",  # Coral
  "#FF7F6A",  # Living Coral
  "#F28C28",  # Mandarin Orange
  "#F4D03F",  # Reef Yellow
  "#36A66D",  # Emerald
  "#8BC34A",  # Lime Reef
  "#B05CC6",  # Orchid
  "#D9487B",  # Raspberry
  "#C62839",  # Crimson
  "#D8C49A",  # Sandstone
  "#6A7885",  # Slate
  "#424242"   # Charcoal
)

reef_muted <- c(
  "#3A5878",  # Navy
  "#5D8BB8",  # Steel Blue
  "#6CB6B6",  # Muted Teal
  "#B8DDD8",  # Pale Aqua
  "#D98C8C",  # Dusty Coral
  "#D79A88",  # Rose Clay
  "#C78B55",  # Terracotta
  "#D8C46D",  # Warm Sand
  "#76A77A",  # Sage Green
  "#9CAD66",  # Olive Green
  "#A58AB8",  # Dusty Lavender
  "#B67C92",  # Dusty Rose
  "#A95D63",  # Brick Red
  "#DDD2B6",  # Light Sand
  "#8A98A6",  # Blue Grey
  "#5B5B5B"   # Graphite
)
# 
# # # Note: we want that small p-values are also small empirical p-values
# # 
# # # function to calculate emprical p-values
# # add_empirical_p <- function(x) {
# #   lapply(x, function(df) {
# #     df$p <- PicMin:::EmpiricalPs(
# #       df$combined_p,
# #       large_i_small_p = FALSE
# #     )
# #     return(df)
# #   })
# # }
# # 
# # # function to create the picmin input r-object list (picmin expect the emp. p-values to be called "p")
# # extract_picmin_inputs <- function(x) {
# #   lapply(x, function(df) {
# #     df %>%
# #       filter(
# #         gene != ""
# #       ) %>%
# #       select(gene, p)
# #   })
# # }
# # 
# # # apply function for all input lists to only contain gene_id and p_combined
# # for (prefix in prefix_fish) {
# #   # define input and output (name)
# #   input_name  <- paste0(prefix, "_gene_results")
# #   output_name <- paste0(prefix, "_picmin_input")
# #   # run functions
# #   result <- get(input_name) |>
# #     add_empirical_p() |>
# #     extract_picmin_inputs()
# #   # define output
# #   assign(output_name, result)
# # }
# # 
# # # combine all separate fish list into an overall list
# # picmin_inputs_all <- setNames(
# #   lapply(prefix_fish, function(prefix) {
# #     get(paste0(prefix, "_picmin_input"))
# #   }),
# #   prefix_fish
# # )
# 
# # old one by oliver
# ## run PicMin for every environmental variable
# for (env in names(PVALS_ALL)) { # for every environmental variable...
#   
#   # set p-values table for current environmental variable
#   all_lins_p = PVALS_ALL[[env]]
#   rownames(all_lins_p) = all_lins_p$gene
#   all_lins_p = all_lins_p[,-1]
#   
#   for (col_i in seq_len(ncol(all_lins_p))) {
#     all_lins_p[, col_i] <- PicMin:::EmpiricalPs(all_lins_p[, col_i], large_i_small_p = FALSE)
#   }
#   
#   ### Run picmin
#   my_picmin_output = RunPicmin(all_lins_p, nullP = nullP)
#   
#   print(env)
#   print(sum(my_picmin_output$pooled_q<0.25, na.rm=T))
#   
#   ### Save results
#   save(my_picmin_output, file=paste0(env,'_my_picmin_output.rda'))
#   
# }
# 
# # record ending time
# t1 = Sys.time()
# # check runtime
# t1-t0
