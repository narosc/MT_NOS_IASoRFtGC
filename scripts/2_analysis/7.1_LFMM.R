#=============================================================================
### seventh analysis 1 - GEA - Latent Factor Mixed Model
#=============================================================================

# This script computes a Genotype Environment Association (GEA) using a Latent 
# Factor Mixed Model (LFMM). 
# Inputs:
# 1) RECIFs_samples -> avg & sd for sample reef cells
# 2) adjusted run table with metadata for samples
# 3) genlight object with imputed/filtered genotype (gt) matrix
# Outputs: 
# 1. LFMM results as r-object list (SNP, z, GIF, p_uncorr, p_gif, qval) 
# 2. List of significantly associated SNPs env. variable as r-object

# Author: Naroa Olivia Schweizer
# last update: 25.01.26

# set correct working directory: setwd("C:/Users/Admin/Documents/Uni_Zuerich/Masters_Thesis/CoralReef_fish_local")

# load libraries
library(ggplot2)        # general-purpose plotting
library(dplyr)          # data wrangling (filter, summarize, join)
library(lfmm)           # run lfmm using this package
library(adegenet)       # PCA & multivariate analysis of genetic data
library(patchwork)      # combining multiple plots

#-----------------------------------------------------------------------------
## define paths and load files
#-----------------------------------------------------------------------------

## sample
# define path to run table
RunTable_path <-"data/AmphiBic_RunTable.csv"

# read run table
RunTable <- read.csv(RunTable_path, header = TRUE, stringsAsFactors = FALSE)

## RECIFs
# define path to filtered recifs data for PCA (avg and sd)
RECIFs_sample_path <- "intermed_outputs/amphiprion_bicinctus/env_data_inputs/AmphiBic_RECIFs_sample.csv"

# read filtered recifs data for PCA (avg and sd data for begining to 2011)
RECIFs_sample <- read.csv(RECIFs_sample_path, header = TRUE)

# add sample row names to env. data sample
rownames(RECIFs_sample) <- RunTable$Run

## Genotype matrix
# define path to filtered and imputed genlight object
genlight_obj_path <- "intermed_outputs/amphiprion_bicinctus/var/amphibic_genlight_filtered_imputed.rds"

# read genlight object
genlight_obj <- readRDS(genlight_obj_path)

## outputs
# create output direcotry for lfmm outputs
dir.create('intermed_outputs/amphiprion_bicinctus/lfmm', showWarnings = FALSE)
# define output path to save etc.
output_path <- "intermed_outputs/amphiprion_bicinctus/lfmm"

#=============================================================================
## 1) prep the input data
#=============================================================================

## genotype matrix
# load genlight object as df
gt <- as.data.frame(as.matrix(genlight_obj))

## environmental matrix
# exclude geographic columns (should be 32 columns)
num_env_data <- RECIFs_sample %>%
  select(-sample_lat, -sample_lon, -allDB_coords.lat, -allDB_coords.lon)

# check if all data (gt and env.) is numeric -> needed for LFMM -> str(df)

#=============================================================================
## 2) calculate LFMM and visualizations
#=============================================================================

# define population structure (how many populations)
K <- 3
# for Petrois volitans -> K = 2
# for Amphiprion bicinctus -> K = 3 (what LEA:snmf least cross entropy)
# for Siphamia tubifer -> K = 1
# for Epinephelus striatus -> K = 1

## two penalty methods for lfmm package
# ridge: shrinks coefficients, keep all predictors, stabilize effective size 
# estimates -> avoid false negatives -> better suited for highly correlated 
# variables, hence polygenic data
# lasso: shrinks and select variables, set some coefficients to zero -> better 
# if only few predictor variables drive variation

#-----------------------------------------------------------------------------
## 2.1) calculate lfmm for each env. variables (sd & avg) separately
#-----------------------------------------------------------------------------

# this part will extract the GIF (genomic inflation factor) for each env
# variable and correct the p-values with the GIF within the looped function.

# create empty list to store GIF results
gif_list <- list()

# number of environmental variables
env_vars <- colnames(num_env_data)

# create list for results of lfmm
results_list <- list()

# loop around each env variable separately
for (var in env_vars) {
  X_var <- as.matrix(num_env_data[[var]])
  
  # fit LFMM
  lfmm_fit <- lfmm_ridge(
    Y = gt,
    X = X_var,
    K = K
  )
  
  # get test results (includes z-scores + GIF)
  # -> calibrate = "gif" so LFMM internally adjusts for GIF
  test_res <- lfmm_test(
    Y = gt,
    X = X_var,
    lfmm = lfmm_fit,
    calibrate = "gif"
  )
  
  # extract z-scores
  zscores <- test_res$score[,1]
  
  # extract GIF (scalar)
  gif_val <- as.numeric(test_res$gif)
  
  # extract uncorrected p-values (raw LFMM p-values)
  p_uncorr <- test_res$pvalue[,1]
  
  # extract GIF-adjusted p-values directly from LFMM
  p_gif <- test_res$calibrated.pvalue[,1]
  
  # save to your results list
  results_list[[var]] <- data.frame(
    SNP = colnames(gt),
    z = zscores,
    GIF = gif_val,
    p_uncorr = p_uncorr,
    p_gif = p_gif
  )
  
  # also store GIF alone for plotting
  gif_list[[var]] <- gif_val
}

#-----------------------------------------------------------------------------
## 2.2) plot GIF and corrected p-value for each variable
#-----------------------------------------------------------------------------

## GIF
# convert to tidy data frame
gif_df <- data.frame(
  envVar = names(gif_list),
  GIF = unlist(gif_list)
)

# plot GIF for each PC
ggplot(gif_df, aes(x = envVar, y = GIF)) +
  geom_col(fill = "grey70") +
  geom_hline(yintercept = 1, linetype = "dashed", color = "black") +
  geom_text(aes(label = round(GIF, 2)), 
            vjust = -0.3, size = 4) +
  
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
  ) +
  
  labs(
    x = "Env Variable",
    y = "Genomic Inflation Factor (GIF)",
    title = "GIF per env. Variable for LFMM"
  )

## p-value
# function to create qq-plots
qq_pvalue_plot <- function(p_uncorr, p_corr, envVar_clean) {
  
  # remove NA
  p_uncorr <- p_uncorr[!is.na(p_uncorr)]
  p_corr   <- p_corr[!is.na(p_corr)]
  
  # sort observed
  p_unc_sorted  <- sort(p_uncorr)
  p_corr_sorted <- sort(p_corr)
  
  # expected p-values
  p_exp <- seq(1/(length(p_unc_sorted)+1), 1, length.out = length(p_unc_sorted))
  
  df <- data.frame(
    expected = -log10(p_exp),
    observed_unc = -log10(p_unc_sorted),
    observed_corr = -log10(p_corr_sorted)
  )
  
  ggplot(df, aes(x = expected)) +
    geom_line(aes(y = expected), linetype = "dashed", color = "black") +
    geom_point(aes(y = observed_unc, color = "Uncorrected"), alpha = 0.6, size = 1.2) +
    geom_point(aes(y = observed_corr, color = "GIF-corrected"), alpha = 0.6, size = 1.2) +
    scale_color_manual(values = c(
      "Uncorrected" = "#A33F39",
      "GIF-corrected" = "#38C3C9"
    )) +
    
    # clean the plot
    theme_minimal(base_size = 11) +
    theme(
      legend.position = "none", 
      plot.title = element_text(size = 10, face = "bold")
    ) +
    
    labs(
      title = envVar_clean,
      x = "Expected",
      y = "Observed"
    )
}

# create list to store all plots in
qq_list <- list()   
# define i = 1
i <- 1

# function to create plot for each env variable
for (var in names(results_list)) {
  
  # remove prefix "s_allDB_"
  var_clean <- gsub("^s_allDB_", "", var)
  
  p_uncorr <- results_list[[var]]$p_uncorr
  p_corr   <- results_list[[var]]$p_gif
  
  qq_list[[i]] <- qq_pvalue_plot(
    p_uncorr = p_uncorr,
    p_corr   = p_corr,
    envVar_clean = var_clean
  )
  
  i <- i + 1
}

# create function legend for patchwork plot
legend_plot <- ggplot(data.frame(x = 1, y = 1)) +
  geom_point(aes(x, y, color = "Uncorrected")) +
  geom_point(aes(x, y, color = "GIF-corrected")) +
  scale_color_manual(values = c("Uncorrected" = "#A33F39",
                                "GIF-corrected" = "#38C3C9")) +
  theme_void() +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(size = 12)
  )

# extract legend for patchwork plot
get_legend <- function(plot) {
  g <- ggplotGrob(plot)
  legend <- g$grobs[which(sapply(g$grobs, function(x) x$name) == "guide-box")][[1]]
  return(legend)
}
# name -> shared legend
shared_legend <- get_legend(legend_plot)

# create patchwork plot 8-pictures in 4 rows, for all the env. variables
qq_patch <- wrap_plots(qq_list, ncol = 8, nrow = 4)

# add top title + legend at bottom
qqplot_comb <- 
  plot_spacer() +
  ggplot() + 
  theme_void() +
  theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5)) +
  labs(title = "QQ-plots for all Env. Variables: p-values uncorrected vs. GIF-corrected in -log10 scale") +
  qq_patch +
  shared_legend +
  plot_layout(
    ncol = 1,
    heights = c(0.1, 0.05, 1, 0.15)
  )

# display
qqplot_comb

#-----------------------------------------------------------------------------
## 2.3) implement correction for multiple testing
#-----------------------------------------------------------------------------

# calculate FDR-adjusted p-values for all variables in result list (BH procedure)
for(var in names(results_list)) {
  results_list[[var]]$qval <- p.adjust(results_list[[var]]$p_gif, method = "BH")
}

# create list to store significantly associated SNPs
sig_snps_list <- list()
# filter for SNPs with a q-value below 0.1 (only keep the most significant)
for(var in names(results_list)) {
  sig_snps_list[[var]] <- results_list[[var]]$SNP[results_list[[var]]$qval < 0.05]
}

## visualize amount of SNPs associated with env. variable
# count associated snps per env. variable and convert to df
sig_snps_counts <- data.frame(
  envVar = names(sig_snps_list),
  n_sig_snps = sapply(sig_snps_list, length)
)

# plot
ggplot(sig_snps_counts, aes(x = envVar, y = n_sig_snps)) +
  geom_col(fill = "grey40") +
  geom_text(aes(label = ifelse(n_sig_snps > 0, n_sig_snps, "")), 
            vjust = -0.3,
            size = 3) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(
    x = "Env Variable",
    y = "# SNPs (q-value < 0.05)",
    title = "# SNPs Sig. Associated with Environmental Variables"
  )

#-----------------------------------------------------------------------------
## 2.4) create manhattan plot
#-----------------------------------------------------------------------------

# create empty list to store Manhattan plots
manhattan_list <- list()

# loop through all environmental variables
for (var in names(results_list)) {
  
  res_df <- results_list[[var]] %>%
    mutate(
      snp_index = 1:n(),                 # sequential SNP index
      logp = -log10(p_gif),              # -log10(GIF-corrected p-value)
      q_signif = qval < 0.1              # highlight significant SNPs
    )
  
  var_clean <- gsub("^s_allDB_", "", var)  # remove prefix for title
  
  manhattan_list[[var]] <- ggplot(res_df, aes(x = snp_index, y = logp, color = q_signif)) +
    geom_point(alpha = 0.7, size = 1.5) +
    scale_color_manual(values = c("FALSE" = "grey70", "TRUE" = "#A33F39")) +
    geom_hline(yintercept = -log10(0.1), linetype = "dashed", color = "black") +
    theme_minimal(base_size = 12) +
    labs(
      x = NULL,  # remove x-axis label for all
      y = NULL,  # remove y-axis label for all
      title = var_clean
    ) +
    theme(
      legend.position = "none",
      axis.text.x = element_text(size = 8),
      axis.text.y = element_text(size = 8)
    )
}

## combine all Manhattan plots into a single figure
# use plot_layout to collect axes, then add a single x/y label with plot_annotation
final_manhattan_plot <- wrap_plots(manhattan_list, ncol = 4) +
  plot_annotation(
    title = "Manhattan Plots for All Environmental Variables",
    theme = theme(plot.title = element_text(size = 16, face = "bold"))
  ) &
  theme(
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12)
  )

# add shared x and y axes
final_manhattan_plot <- final_manhattan_plot &
  plot_layout(guides = "collect")  

# display
final_manhattan_plot

#=============================================================================
## 3) save important data
#=============================================================================

## LFMM results
# save LFMM results as r-object (SNP, z, GIF, p_uncorr, p_gif, qval)
save(results_list, file = "intermed_outputs/amphiprion_bicinctus/lfmm/lfmm_results.RData")

# save list of sig. associated SNPs to env. variable as r-object
save(sig_snps_list, file = "intermed_outputs/amphiprion_bicinctus/lfmm/lfmm_sigAsso_SNPsEnvVar.RData")
