#=============================================================================
### seventh analysis 2 - GEA - Correlation Analysis
#=============================================================================

# This script computes a Genotype Environment Association (GEA) using a plain
# correlation analysis (CA).
# Inputs:
# 1) RECIFs_samples -> avg & sd for sample reef cells
# 2) adjusted run table with metadata for samples
# 3) genlight object with imputed/filtered genotype (gt) matrix
# Outputs: 
# 1. CA results as r-object list (SNP, pvalue, tau, qval) 
# 2. List of significantly associated SNPs env. variable as r-object

# Author: Naroa Olivia Schweizer
# last update: 10.06.26

# set correct working directory: setwd("C:/Users/Admin/Documents/Uni_Zuerich/Masters_Thesis/CoralReef_fish_local")

# load libraries
library(ggplot2)        # general-purpose plotting
library(dplyr)          # data wrangling (filter, summarize, join)
library(adegenet)       # PCA & multivariate analysis of genetic data
library(patchwork)      # combining multiple plots

#-----------------------------------------------------------------------------
## define paths and load files
#-----------------------------------------------------------------------------

## sample
# define path to run table
RunTable_path <-"data/EpiStri_RunTable.csv"

# read run table
RunTable <- read.csv(RunTable_path, header = TRUE, stringsAsFactors = FALSE)

## RECIFs
# define path to filtered recifs data for PCA (avg and sd)
RECIFs_sample_path <- "intermed_outputs/epinephelus_striatus/env_data_inputs/EpiStri_RECIFs_sample.csv"

# read filtered recifs data for PCA (avg and sd data for begining to 2011)
RECIFs_sample <- read.csv(RECIFs_sample_path, header = TRUE)
# add sample row names to env. data sample
rownames(RECIFs_sample) <- RunTable$Run

## Genotype matrix
# define path to filtered and imputed genlight object
genlight_obj_path <- "intermed_outputs/epinephelus_striatus/var/epistri_genlight_filtered_imputed.rds"

# read genlight object
genlight_obj <- readRDS(genlight_obj_path)

## outputs
# create output direcotry for lfmm outputs
dir.create('intermed_outputs/epinephelus_striatus/ca', showWarnings = FALSE)
# define output path to save etc.
output_path <- "intermed_outputs/epinephelus_striatus/ca"

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

#=============================================================================
## 2) calculate correlation analysis
#=============================================================================

# create empty list
results_list <- list()

# loop over environmental variables
for (env_name in colnames(num_env_data)) {
  # extraxt env variable
  env <- num_env_data[[env_name]]
  # create result object
  env_results <- vector(
    "list",
    ncol(gt)
  )
  
  counter <- 1
  # calculate correlation for each SNP with each environment
  for (snp_name in colnames(gt)) {
    
    snp <- gt[, snp_name]
    # kendall: measures the strength of monotonic association between two sets of ranked or ordinal data
    test <- suppressWarnings(
      cor.test(
        snp,
        env,
        method = "kendall"
      )
    )
    # paste SNP name, tau (no corr = 0, 1 to -1 correlation), p-value
    env_results[[counter]] <- data.frame(
      SNP = snp_name,
      tau = unname(test$estimate),
      pvalue = test$p.value
    )
    
    counter <- counter + 1
  }
  # create result list 
  results_list[[env_name]] <- do.call(
    rbind,
    env_results
  )
}

#-----------------------------------------------------------------------------
## 2.1) implement FDR, calculate q-values
#-----------------------------------------------------------------------------

# calculate q-value (Benjamin-Hochberg)
for(var in names(results_list)) {
  
  results_list[[var]]$qval <- p.adjust(
    results_list[[var]]$pvalue,
    method = "BH"
  )
  
}

#=============================================================================
## 3) create visualisations
#=============================================================================

#-----------------------------------------------------------------------------
## 3.1) qq-plots
#-----------------------------------------------------------------------------

# function to create qq plots
qq_pvalue_plot <- function(pvals, envVar_clean) {
  pvals <- pvals[!is.na(pvals)]
  p_obs <- sort(pvals)
  p_exp <- seq(
    1/(length(p_obs)+1),
    1,
    length.out = length(p_obs)
  )
  df <- data.frame(
    expected = -log10(p_exp),
    observed = -log10(p_obs)
  )
  ggplot(df, aes(x = expected, y = observed)) +
    geom_abline(
      slope = 1,
      intercept = 0,
      linetype = "dashed",
      color = "grey30",
      linewidth = 0.8
    ) +
    geom_point(
      color = "#A33F39",
      alpha = 0.6,
      size = 1
    ) +
    theme_minimal(base_size = 11) +
    labs(
      title = envVar_clean,
      x = "Expected -log10(p)",
      y = "Observed -log10(p)"
    )
}

# create list with plots to store
qq_list <- list()
# define i
i <- 1
# loop over all variables
for(var in names(results_list)) {
  qq_list[[i]] <- qq_pvalue_plot(
    pvals = results_list[[var]]$pvalue,
    envVar_clean = var
  )
  i <- i + 1
}

# combine plots
qq_patch <- wrap_plots(
  qq_list,
  ncol = 8
)
# plot 
qq_patch +
  plot_annotation(
    title = "QQ-plots for all environmental variables Correlation Analysis",
    theme = theme(plot.title = element_text(size = 16, face = "bold"))
  ) &
  theme(
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12)
  )

#-----------------------------------------------------------------------------
## 3.2) manhattan plots
#-----------------------------------------------------------------------------

# create function to calculate manhattan plots
manhattan_list <- list()
# run loop
for(var in names(results_list)) {
  res_df <- results_list[[var]] %>%
    mutate(
      snp_index = 1:n(),
      logp = -log10(pvalue),
      q_signif = qval < 0.25
    )
  manhattan_list[[var]] <-
    ggplot(
      res_df,
      aes(
        x = snp_index,
        y = logp,
        color = q_signif
      )
    ) +
    geom_point(
      alpha = 0.7,
      size = 1.3
    ) +
    scale_color_manual(
      values = c(
        "FALSE" = "grey70",
        "TRUE" = "#A33F39"
      )
    ) +
    geom_hline(
      yintercept = -log10(0.05),
      linetype = "dashed"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      legend.position = "none"
    ) +
    labs(
      title = var,
      x = NULL,
      y = NULL
    )
}

# create plots
final_manhattan_plot <-  wrap_plots(manhattan_list, ncol = 4
  ) +
  plot_annotation(
    title = "Manhattan Plots for all Environmental Variables Correlation Analysis",
    theme = theme(plot.title = element_text(size = 16, face = "bold"))
  ) &
  theme(
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12)
  )
# plot
final_manhattan_plot

#-----------------------------------------------------------------------------
## 3.3) tau p-value distribution plot
#-----------------------------------------------------------------------------

# create list for tau plots
tau_sig_list <- list()
# loop over all env variables
for(var in names(results_list)) {
  var_clean <- gsub(
    "^s_allDB_",
    "",
    var
  )
  tau_sig_list[[var]] <-
    ggplot(
      results_list[[var]],
      aes(
        x = tau,
        y = -log10(pvalue),
        color = qval < 0.25
      )
    ) +
    geom_point(
      alpha = 0.5,
      size = 1.2
    ) +
    scale_color_manual(
      values = c(
        "FALSE" = "grey40",
        "TRUE" = "#A33F39"
      )
    ) +
    geom_vline(
      xintercept = 0,
      linetype = "dashed"
    ) +
    geom_hline(
      yintercept = -log10(0.01),
      linetype = "dashed"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      legend.position = "none",
      plot.title = element_text(
        size = 10,
        face = "bold"
      )
    ) +
    labs(
      title = var_clean,
      x = NULL,
      y = NULL
    )
}

# create patchwork plot
tau_sig_patch <- wrap_plots(
  tau_sig_list,
  ncol = 8,
  nrow = 4
) +
  plot_annotation(
    title = "Tau-P-Value Distribution Plots for all Environmental Variables Correlation Analysis",
    theme = theme(plot.title = element_text(size = 16, face = "bold"))
  ) &
  theme(
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12)
  )

# display plot
tau_sig_patch

#=============================================================================
## 3) save important data
#=============================================================================

## CA results
# save LFMM results as r-object (SNP, tau, pvalue, qval)
save(results_list, file = "intermed_outputs/epinephelus_striatus/ca/ca_results.RData")
