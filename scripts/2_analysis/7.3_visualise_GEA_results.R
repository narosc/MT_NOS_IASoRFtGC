#=============================================================================
### 7.3 - visulaise GEA results - CA and LFMM - overall figure
#=============================================================================

# This script creates a summary figure showing the GEA (LFMM and CA) results.
# Input:
# 1) for all 4 species: lfmm_results.RData and ca_results.RData
# Output: 
# 1. summary figure visualising GEA results

# set correct working directory: 
setwd("C:/Users/Admin/Documents/Uni_Zuerich/Masters_Thesis/CoralReef_fish_local")

# Author: Naroa Olivia Schweizer
# last update: 16.08.2026

# load libraries
library(dplyr)           # data wrangling (filter, summarize, join)
library(ggplot2)         # general-purpose plotting
library(tidyr)           # reshape and tidy data frames

#-----------------------------------------------------------------------------
## define paths and load files
#-----------------------------------------------------------------------------

## Pterois volitans
# define path to r-object LFMM
SNPs_list_path <- "intermed_outputs/pterois_volitans/lfmm/lfmm_results.RData"

# load r-object (results_list)
load(SNPs_list_path)
ptevol_lfmm <- results_list

# define path to r-object CA
SNPs_list_path <- "intermed_outputs/pterois_volitans/ca/ca_results.RData"

# load r-object (results_list)
load(SNPs_list_path)
ptevol_ca <- results_list

## Amphiprion bicinctus
# define path to r-object LFMM
SNPs_list_path <- "intermed_outputs/amphiprion_bicinctus/lfmm/lfmm_results.RData"

# load r-object (results_list)
load(SNPs_list_path)
amphibic_lfmm <- results_list

# define path to r-object CA
SNPs_list_path <- "intermed_outputs/amphiprion_bicinctus/ca/ca_results.RData"

# load r-object (results_list)
load(SNPs_list_path)
amphibic_ca <- results_list

## Siphamia tubifer
# define path to r-object LFMM
SNPs_list_path <- "intermed_outputs/siphamia_tubifer/lfmm/lfmm_results.RData"

# load r-object (results_list)
load(SNPs_list_path)
siphtub_lfmm <- results_list

# define path to r-object CA
SNPs_list_path <- "intermed_outputs/siphamia_tubifer/ca/ca_results.RData"

# load r-object (results_list)
load(SNPs_list_path)
siphtub_ca <- results_list

## Epinephelus striatus
# define path to r-object LFMM
SNPs_list_path <- "intermed_outputs/epinephelus_striatus/lfmm/lfmm_results.RData"

# load r-object (results_list)
load(SNPs_list_path)
epistri_lfmm <- results_list

# define path to r-object CA
SNPs_list_path <- "intermed_outputs/epinephelus_striatus/ca/ca_results.RData"

# load r-object (results_list)
load(SNPs_list_path)
epistri_ca <- results_list

# create vectors containing prefix names of all species
prefix_fish <- c("ptevol","amphibic", "siphtub", "epistri")

# define output path for the summary stats csv (adjust as needed)
summary_save_path <- "outputs/GEA_visualisation/GEA_summary_stats.csv"

# define output path for the summary stats csv (adjust as needed)
summary_save_path2 <- "outputs/GEA_visualisation/GEA_total_counts_perSandEnv.csv"

#=============================================================================
## 1) create concordance bar plot for visualisation
#=============================================================================

#-----------------------------------------------------------------------------
## 1.1) specifics only for Epinephelus striatus
#-----------------------------------------------------------------------------

# SNPs to exclude for Epinephelus striatus (heterozygous but monomorphic -
# survived filtering but shouldn't be included)
NA_NaN_snps <- c(
  "CM069294.1_12128051",
  "CM069295.1_2820193",
  "CM069296.1_2928336",
  "CM069297.1_4068134",
  "CM069298.1_40857250",
  "CM069299.1_7010888",
  "CM069300.1_33637199",
  "CM069301.1_12515129",
  "CM069302.1_42396464",
  "CM069303.1_22455541",
  "CM069304.1_106062",
  "CM069305.1_5136027",
  "CM069306.1_12026293",
  "CM069307.1_3857487",
  "CM069308.1_925475",
  "CM069309.1_35744269",
  "CM069310.1_33416869",
  "CM069311.1_23865501",
  "CM069312.1_32003298",
  "CM069313.1_27026024",
  "CM069316.1_52950",
  "CM069317.1_19741987"
)

# function to remove these SNPs from every env-variable data frame within
# a species' result list (works for both lfmm and ca structures, since
# both have a "SNP" column)
remove_snps <- function(result_list, snps_to_remove) {
  lapply(result_list, function(df) {
    df[!(df$SNP %in% snps_to_remove), ]
  })
}

# apply to epistri only
epistri_lfmm <- remove_snps(epistri_lfmm, NA_NaN_snps)
epistri_ca   <- remove_snps(epistri_ca, NA_NaN_snps)

#-----------------------------------------------------------------------------
## 1.2) prepare plot df (long format)
#-----------------------------------------------------------------------------

# create long formate df for plotting
build_snp_sig_table <- function(lfmm_list, ca_list, threshold = q_threshold) {
  lfmm_long <- bind_rows(lapply(names(lfmm_list), function(env) {
    df <- lfmm_list[[env]]
    data.frame(SNP = df$SNP, env_var = env, sig_lfmm = df$qval < threshold)
  }))
  
  ca_long <- bind_rows(lapply(names(ca_list), function(env) {
    df <- ca_list[[env]]
    data.frame(SNP = df$SNP, env_var = env, sig_ca = df$qval < threshold)
  }))
  
  full_join(lfmm_long, ca_long, by = c("SNP", "env_var"))
}

# create species list
species_lists <- list(
  ptevol   = list(lfmm = ptevol_lfmm,   ca = ptevol_ca),
  amphibic = list(lfmm = amphibic_lfmm, ca = amphibic_ca),
  siphtub  = list(lfmm = siphtub_lfmm,  ca = siphtub_ca),
  epistri  = list(lfmm = epistri_lfmm,  ca = epistri_ca)
)

# define q-value threshold
q_threshold <- 0.05

# define species lables
species_labels <- c(
  ptevol   = "Pterois volitans",
  amphibic = "Amphiprion bicinctus",
  siphtub  = "Siphamia tubifer",
  epistri  = "Epinephelus striatus"
)
# run function
all_sig_tables <- lapply(names(species_lists), function(sp) {
  build_snp_sig_table(species_lists[[sp]]$lfmm, species_lists[[sp]]$ca) %>%
    mutate(species = sp)
})
# clean env. variable names (remove prefix)
all_sig_tables <- bind_rows(all_sig_tables) %>%
  mutate(env_var_clean = sub("^s_allDB_", "", env_var))

# create per-method SNP counts
concordance_df <- all_sig_tables %>%
  group_by(species, env_var_clean) %>%
  summarise(
    LFMM = sum(sig_lfmm, na.rm = TRUE),
    CA   = sum(sig_ca, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_longer(c(LFMM, CA), names_to = "category", values_to = "n_snps") %>%
  mutate(
    category = factor(category, levels = c("LFMM", "CA")),
    species = recode(species, !!!species_labels)
  )

#-----------------------------------------------------------------------------
## 1.3) create concordance summary figure
#-----------------------------------------------------------------------------

# define fill colurs
fill_cols <- c(
  LFMM = "#6CB6B6",
  CA   = "#F46D75"
)
# define outline colours (slightly darker) - also used for the SNP count labels
border_cols <- c(
  LFMM = "#408280",
  CA   = "#B20E16"
)

# shared dodge width, so bars and their count labels line up
dodge_width <- 0.9

# force axis breaks to whole numbers only
integer_breaks <- function(n = 5) {
  function(x) unique(round(pretty(x, n)))
}

# create summary plot
p_concordance <- ggplot(
  concordance_df,
  aes(
    x = env_var_clean,
    y = n_snps,
    fill = category,
    colour = category
  )
) +
  geom_col(
    position = position_dodge(width = dodge_width),
    width = 0.95,   # thicker bars (each rendered bar = width / n_groups)
    linewidth = 0.3
  ) +
  geom_text(
    aes(label = ifelse(n_snps >= 1, n_snps, "")),   # only label bars with >= 1 SNP
    position = position_dodge(width = dodge_width),
    hjust = -0.2,
    size = 2.7,
    show.legend = FALSE
  ) +
  scale_fill_manual(
    values = fill_cols,
    name = "GEA Method",
    labels = c(
      LFMM = "LFMM",
      CA   = "CA"
    ),
    guide = guide_legend(title.position = "top",
                         title.hjust = 0.5)
  ) +
  scale_colour_manual(
    values = border_cols,
    guide = "none"   # no second legend - colour is only used for the outline/labels
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.15)),   # room for the count labels
    breaks = integer_breaks(),
    labels = scales::label_number(accuracy = 1)
  ) +
  coord_flip() +
  facet_wrap(~species, scales = "free_x") +
  labs(
    x = NULL,
    y = paste0("SNP count (q < ", q_threshold, ")"),
    title = "GEA Summary Figure"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    strip.text = element_text(face = "bold",
                              size = 12),
    legend.position = "bottom",
    legend.title = element_text(
      face = "bold",
      hjust = 0.5,
      size = 10
    ),
    plot.title = element_text(
      hjust = 0.5,
      face = "bold",
      size = 18
    ),
    axis.text.y = element_text(size = 8),
    axis.text.x = element_text(size = 8)
  )

# display
p_concordance
# saved 750 x 1000

#=============================================================================
## 2) summary statistics for thesis text (median / IQR)
#=============================================================================

#-----------------------------------------------------------------------------
## 2.1) median/IQR SNP count per species and GEA method
#-----------------------------------------------------------------------------

# calculate median/IQR SNP count per environmental variable
# separately for each species and GEA method
summary_stats <- concordance_df %>%
  group_by(species, category, env_var_clean) %>%
  summarise(
    n_snps = sum(n_snps),
    .groups = "drop"
  ) %>%
  group_by(species, category) %>%
  summarise(
    median_n_snps = round(median(n_snps), 1),
    IQR_n_snps    = round(IQR(n_snps), 1),
    .groups = "drop"
  ) %>%
  rename(
    group = species,
    method = category
  ) %>%
  mutate(
    grouping = "species"
  ) %>%
  select(grouping, group, method, median_n_snps, IQR_n_snps) %>%
  arrange(group, method)

#-----------------------------------------------------------------------------
## 2.2) overall method comparison
#-----------------------------------------------------------------------------

# median/IQR SNP count per method and environmental variable,
# summed across species
method_summary <- concordance_df %>%
  group_by(category, env_var_clean) %>%
  summarise(
    n_snps_total = sum(n_snps),
    .groups = "drop"
  ) %>%
  group_by(category) %>%
  summarise(
    median_n_snps = round(median(n_snps_total), 1),
    IQR_n_snps    = round(IQR(n_snps_total), 1),
    .groups = "drop"
  ) %>%
  rename(
    group = category
  ) %>%
  mutate(
    grouping = "method",
    method = group
  ) %>%
  select(grouping, group, method, median_n_snps, IQR_n_snps)

#-----------------------------------------------------------------------------
## 2.3) combine and save
#-----------------------------------------------------------------------------

# combine species-specific and overall method summaries
summary_stats <- bind_rows(
  method_summary,
  summary_stats
)

# inspect
summary_stats

# save
write.csv(summary_stats, file = summary_save_path, row.names = FALSE)

# save total SNPs per env. and species
write.csv(concordance_df, file = summary_save_path2, row.names = FALSE)
