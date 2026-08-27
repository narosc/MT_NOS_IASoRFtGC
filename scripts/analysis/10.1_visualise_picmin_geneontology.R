#==============================================================================
### tenth - visualisation of picmin and gene ontology
#==============================================================================

# This script visualises the picmin outputs by also combining them with the
# gene ontology outputs.
# Inputs:
# 1) picmin output -> results list (CA and LFMM)
# 2) gene ontology outputs -> BP/CC/MF_gsc for CA and LFMM
# Output:
# 1. main summary plot: panel a + panel b (picmin hits), patchwork side by side
# 2. supplementary plot: panel c (GO heatmap), standalone

# set correct working directory:
setwd("C:/Users/Admin/Documents/Uni_Zuerich/Masters_Thesis/CoralReef_fish_local")

# Author: Naroa Olivia Schweizer
# last update: 06.08.26

# load libraries
library(dplyr)       # data wrangling (filter, summarise, join)
library(tidyr)       # reshape/tidy data frames
library(ggplot2)     # plotting
library(patchwork)   # combine bar charts + heatmap into one figure
library(ggnewscale)  # lets LFMM and CA each keep their own fill scale/legend within one panel
library(SetRank)     # GO enrichment test, built for buildSetCollection() gsc objects
library(igraph)      # to pull the results table out of the SetRank network object

#-----------------------------------------------------------------------------
## define paths and load files
#-----------------------------------------------------------------------------

## PicMin outputs
# lfmm
picmin_path <- "intermed_outputs/picmin/picmin_results_lfmm.rda"
load(picmin_path)                       # object: picmin_results
lfmm_picmin <- picmin_results

# ca
picmin_path <- "intermed_outputs/picmin/picmin_results_ca.rda"
load(picmin_path)                       # object: picmin_results
ca_picmin <- picmin_results

## Gene Ontology output
# biological process
GO_BP_lfmm_path <- "intermed_outputs/setrank/BP_gsc_lfmm.rda"
load(GO_BP_lfmm_path)                   # object: BP_gsc
GO_BP_collection <- BP_gsc              # reused for both methods

# cellular component
GO_CC_lfmm_path <- "intermed_outputs/setrank/CC_gsc_lfmm.rda"
load(GO_CC_lfmm_path)                   # object: CC_gsc
GO_CC_collection <- CC_gsc              # reused for both methods

# molecular function
GO_MF_lfmm_path <- "intermed_outputs/setrank/MF_gsc_lfmm.rda"
load(GO_MF_lfmm_path)                   # object: MF_gsc
GO_MF_collection <- MF_gsc              # reused for both methods

#==============================================================================
## 1) prepare PicMin/GO output for plotting
#==============================================================================

#-----------------------------------------------------------------------------
## 1.1) prepare PicMin lists
#-----------------------------------------------------------------------------

# function to retrieve env variables, selects genes with q-value < 0.25, retruns df in long formate
# -> n_est is NOT filtered here, all contributing-species counts are kept as is
collapse_picmin_list <- function(picmin_list, method_label, q_cutoff = 0.25) {
  rows <- Map(function(df, env_name) {
    df <- df[df$pooled_q < q_cutoff, , drop = FALSE]
    if (nrow(df) == 0) return(NULL)
    # return df with env. variable, method lable, gene id, n_est (contributing species), pooled_q value
    data.frame(
      env_var  = env_name,
      method   = method_label,
      locus    = df$locus,
      n_est    = df$n_est,
      pooled_q = df$pooled_q,
      stringsAsFactors = FALSE
    )
  }, picmin_list, names(picmin_list))
  dplyr::bind_rows(rows)
}

# run function for lfmm and ca
hits_lfmm <- collapse_picmin_list(lfmm_picmin, "LFMM")
hits_ca   <- collapse_picmin_list(ca_picmin,   "CA")
# collaps picmin outputs to one df
hits_all  <- dplyr::bind_rows(hits_lfmm, hits_ca)

## fix the method factor level order once -> LFMM top, CA bottom (for visualisation later)
hits_all$method <- factor(hits_all$method, levels = c("CA", "LFMM"))

#-----------------------------------------------------------------------------
## 1.2) prepare y-axis, legends, titles for plots
#-----------------------------------------------------------------------------

# function to count hits per env. variable, order descending,
# remove env. variables with no hits -> fix for all plots
all_env_vars <- union(names(lfmm_picmin), names(ca_picmin))

env_order <- hits_all %>%
  dplyr::count(env_var, name = "n_total") %>%
  dplyr::right_join(data.frame(env_var = all_env_vars), by = "env_var") %>%
  dplyr::mutate(n_total = tidyr::replace_na(n_total, 0)) %>%
  dplyr::filter(n_total > 0) %>%
  dplyr::arrange(dplyr::desc(n_total)) %>%
  dplyr::pull(env_var)

# enforce this env. order in hits_all df
hits_all$env_var <- factor(hits_all$env_var, levels = env_order)

# get ridd of env_prefix for the variables
strip_env_prefix <- function(x) sub("^s_allDB_?", "", x)

# create consistend legend text/title sizes for all panels
# row_grid_theme <- theme(
#   panel.spacing.y   = unit(0, "lines"),
#   panel.border      = element_rect(colour = "grey65", fill = NA, linewidth = 0.3),
#   panel.grid.major.y = element_blank(),
#   panel.grid.minor   = element_blank(),
#   axis.text.y  = element_blank(),
#   axis.ticks.y = element_blank(),
#   legend.title = element_text(size = 9),
#   legend.text  = element_text(size = 8)
# )
row_grid_theme <- theme(
  panel.spacing.y   = unit(0, "lines"),
  panel.border      = element_rect(colour = "grey65", fill = NA, linewidth = 0.3),
  panel.grid.major.y = element_blank(),
  panel.grid.minor   = element_blank(),
  axis.text.y  = element_blank(),
  axis.ticks.y = element_blank()
)

# center bottom legends
# -> explicit order = keeps CA above LFMM in the legend box, regardless of layer draw order
ca_legend_guide   <- guide_legend(title.position = "top", title.hjust = 0.5, nrow = 1, order = 1)
lfmm_legend_guide <- guide_legend(title.position = "top", title.hjust = 0.5, nrow = 1, order = 2)

#-----------------------------------------------------------------------------
## 1.3) create colour ramps for methods
#-----------------------------------------------------------------------------

# define main lfmm and ca colour
lfmm_colour <- "#6CB6B6"   # turquoise
ca_colour   <- "#F46D75"   # coral

# build sequential shated of one base hue
make_shades <- function(base_colour, n) {
  pale <- grDevices::colorRampPalette(c("white", base_colour))(10)[4]
  deep <- grDevices::colorRampPalette(c(base_colour, "black"))(10)[4]
  grDevices::colorRampPalette(c(pale, base_colour, deep))(n)
}

# darkest tone method's ramp, reused as the colour for the gene count lables
label_colour_lfmm <- make_shades(lfmm_colour, 10)[10]
label_colour_ca   <- make_shades(ca_colour,   10)[10]

#-----------------------------------------------------------------------------
## 1.4) prepare data for genes hit per env with coloured N-species contribution
#-----------------------------------------------------------------------------

# retrieve levels of contributing species (2,3 & 4)
n_est_levels <- sort(unique(hits_all$n_est))
# count genes per method
panelA_counts <- hits_all %>%
  dplyr::count(env_var, method, n_est, name = "n_genes")
# enforce env_order on panel
panelA_grid <- tidyr::expand_grid(
  env_var = factor(env_order, levels = env_order),
  method  = factor(c("CA", "LFMM"), levels = c("CA", "LFMM")),
  n_est   = n_est_levels
)
# join pnale counts?
panelA_data <- panelA_grid %>%
  dplyr::left_join(panelA_counts, by = c("env_var", "method", "n_est")) %>%
  dplyr::mutate(n_genes = tidyr::replace_na(n_genes, 0))

#-----------------------------------------------------------------------------
## 1.5) prepare data for genes hit per env with coloured FDR
#-----------------------------------------------------------------------------

# define q-breaks and labels (exclusive intervals)
# -> rescaled to the 0.25 cutoff, same 5-bin granularity as before (0.05-wide bins)
q_breaks <- c(-Inf, 0.05, 0.1, 0.15, 0.2, 0.25)
q_labels <- c("<0.05", "0.05-0.1", "0.1-0.15", "0.15-0.2", "0.2-0.25")

hits_all <- hits_all %>%
  dplyr::mutate(q_bin = cut(pooled_q, breaks = q_breaks, labels = q_labels, right = FALSE))

# export full hits table (env_var, method, locus, n_est, pooled_q, q_bin) as csv, for the thesis results text
# dir.create("outputs", showWarnings = FALSE)
# hits_all_export <- hits_all %>%
#   dplyr::arrange(env_var, method, pooled_q)
# write.csv(hits_all_export, "outputs/picmin_4s/picmin_hits_summary.csv", row.names = FALSE)

panelB_counts <- hits_all %>%
  dplyr::count(env_var, method, q_bin, name = "n_genes")

panelB_grid <- tidyr::expand_grid(
  env_var = factor(env_order, levels = env_order),
  method  = factor(c("CA", "LFMM"), levels = c("CA", "LFMM")),
  q_bin   = factor(q_labels, levels = q_labels)
)

panelB_data <- panelB_grid %>%
  dplyr::left_join(panelB_counts, by = c("env_var", "method", "q_bin")) %>%
  dplyr::mutate(n_genes = tidyr::replace_na(n_genes, 0))

# create a total bar, same bar for panel a and b -> just colouring different
bar_totals <- panelA_data %>%
  dplyr::group_by(env_var, method) %>%
  dplyr::summarise(total = sum(n_genes), .groups = "drop")

## shared x-axis limit for panels a AND b (same underlying gene totals)
x_max <- max(bar_totals$total)

#-----------------------------------------------------------------------------
## 1.6) prepare data for GO plot
#-----------------------------------------------------------------------------

# set cores (serial here not parallel)
options(mc.cores = 1)

# run set rank
run_go_enrichment_for_variable <- function(env_name, hits_all, go_collection,
                                           set_p_cutoff = 0.01, fdr_cutoff = 0.05) {
  genes <- hits_all %>%
    dplyr::filter(env_var == env_name) %>%
    dplyr::pull(locus) %>%
    unique()
  
  if (length(genes) == 0) return(NULL)
  
  net <- SetRank::setRankAnalysis(genes, go_collection, use.ranks = FALSE,
                                  setPCutoff = set_p_cutoff, fdrCutoff = fdr_cutoff)
  
  if (igraph::vcount(net) == 0) return(NULL)
  
  res <- igraph::as_data_frame(net, what = "vertices")
  data.frame(env_var = env_name,
             go_term = res$description,
             fdr     = res$adjustedPValue,
             stringsAsFactors = FALSE)
}

# run enrichment separately per GO domain (BP/CC/MF), tag rows with the domain they came from
panelC_BP <- dplyr::bind_rows(
  lapply(env_order, run_go_enrichment_for_variable,
         hits_all = hits_all, go_collection = GO_BP_collection)
)
if (nrow(panelC_BP) > 0) panelC_BP$ontology <- "BP"

panelC_CC <- dplyr::bind_rows(
  lapply(env_order, run_go_enrichment_for_variable,
         hits_all = hits_all, go_collection = GO_CC_collection)
)
if (nrow(panelC_CC) > 0) panelC_CC$ontology <- "CC"

panelC_MF <- dplyr::bind_rows(
  lapply(env_order, run_go_enrichment_for_variable,
         hits_all = hits_all, go_collection = GO_MF_collection)
)
if (nrow(panelC_MF) > 0) panelC_MF$ontology <- "MF"

# combine all three GO domains into one df for panel c
panelC_data <- dplyr::bind_rows(panelC_BP, panelC_CC, panelC_MF)

panelC_data$env_var <- factor(panelC_data$env_var, levels = env_order)

# keep ontology in a fixed reading order: BP, CC, MF
panelC_data$ontology <- factor(panelC_data$ontology, levels = c("BP", "CC", "MF"))

# keep GO terms ordered alphabetically along x, within each ontology
panelC_data$go_term <- factor(panelC_data$go_term, levels = sort(unique(panelC_data$go_term)))

#==============================================================================
## 2) create plots and combine into main + supplementary summary plots
#==============================================================================

#-----------------------------------------------------------------------------
## 2.1) plot panel left: gene hits, coloured n_est
#-----------------------------------------------------------------------------

# create lfmm and ca colour gradients
shades_a_lfmm <- setNames(make_shades(lfmm_colour, length(n_est_levels)), n_est_levels)
shades_a_ca   <- setNames(make_shades(ca_colour,   length(n_est_levels)), n_est_levels)

# plot panel left
# -> CA layer/scale added first so CA's legend sits on top (LFMM legend below)
panel_a_plot <-
  ggplot() +
  geom_col(data = dplyr::filter(panelA_data, method == "CA"),
           aes(x = n_genes, y = method, fill = factor(n_est)),
           position = "stack", width = 0.8) +
  scale_fill_manual(values = shades_a_ca, name = "Number of contributing species (CA)", guide = ca_legend_guide) +
  ggnewscale::new_scale_fill() +
  geom_col(data = dplyr::filter(panelA_data, method == "LFMM"),
           aes(x = n_genes, y = method, fill = factor(n_est)),
           position = "stack", width = 0.8) +
  scale_fill_manual(values = shades_a_lfmm, name = "Number of contributing species (LFMM)", guide = lfmm_legend_guide) +
  # total-count labels printed just past the end of each bar
  geom_text(data = dplyr::filter(bar_totals, method == "CA"),
            aes(x = total, y = method, label = total),
            hjust = -0.15, size = 4.5, colour = label_colour_ca) +
  geom_text(data = dplyr::filter(bar_totals, method == "LFMM"),
            aes(x = total, y = method, label = total),
            hjust = -0.15, size = 4.5, colour = label_colour_lfmm) +
  facet_grid(rows = vars(env_var), switch = "y", drop = FALSE,
             labeller = labeller(env_var = strip_env_prefix)) +
  scale_x_continuous(limits = c(0, x_max), expand = expansion(mult = c(0, 0.15))) +
  labs(x = "Number of PicMin Genes (q-value < 0.25)",
       y = NULL,
       title = "Species Contributions to PicMin Genes") +
  theme_minimal(base_size = 11) +
  theme(strip.placement = "outside",
        # bigger env. variable row labels + bigger x-axis title
        strip.text.y.left = element_text(angle = 0, hjust = 1, size = 11),
        axis.title.x = element_text(size = 13),
        axis.text.x = element_text(size = 11),
        legend.title = element_text(size = 11),
        legend.text = element_text(size = 10),
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.box = "vertical",
        legend.box.just = "center",
        legend.box.spacing = unit(8, "pt"),
        legend.spacing.y = unit(1, "pt"),
        legend.margin = margin(t = 15, b = 0),
        plot.margin = margin(t = 5.5, r = 15, b = 5.5, l = 5.5)) +
  row_grid_theme

#-----------------------------------------------------------------------------
## 2.2) plot panel middle: gene hits, coloured FDR
#-----------------------------------------------------------------------------

# create colour gradient CA and LFMM
shades_b_lfmm <- setNames(rev(make_shades(lfmm_colour, length(q_labels))), q_labels)
shades_b_ca   <- setNames(rev(make_shades(ca_colour,   length(q_labels))), q_labels)

# create plot middle
# -> CA layer/scale added first so CA's legend sits on top (LFMM legend below)
panel_b_plot <-
  ggplot() +
  # reverse = TRUE -> lowest q-value bin stacks first, ending up on the left (brighter/higher q to the right)
  geom_col(data = dplyr::filter(panelB_data, method == "CA"),
           aes(x = n_genes, y = method, fill = q_bin),
           position = position_stack(reverse = TRUE), width = 0.8) +
  scale_fill_manual(values = shades_b_ca, name = "Pooled q-value (CA)", guide = ca_legend_guide) +
  ggnewscale::new_scale_fill() +
  geom_col(data = dplyr::filter(panelB_data, method == "LFMM"),
           aes(x = n_genes, y = method, fill = q_bin),
           position = position_stack(reverse = TRUE), width = 0.8) +
  scale_fill_manual(values = shades_b_lfmm, name = "Pooled q-value (LFMM)", guide = lfmm_legend_guide) +
  # total-count labels printed just past the end of each bar
  geom_text(data = dplyr::filter(bar_totals, method == "CA"),
            aes(x = total, y = method, label = total),
            hjust = -0.15, size = 4.5, colour = label_colour_ca) +
  geom_text(data = dplyr::filter(bar_totals, method == "LFMM"),
            aes(x = total, y = method, label = total),
            hjust = -0.15, size = 4.5, colour = label_colour_lfmm) +
  facet_grid(rows = vars(env_var), switch = "y", drop = FALSE,
             labeller = labeller(env_var = strip_env_prefix)) +
  scale_x_continuous(limits = c(0, x_max), expand = expansion(mult = c(0, 0.15))) +
  labs(x = "Number of PicMin Genes (q-value < 0.25)",
       y = NULL,
       title = "Quality Assessment for PicMin Genes") +
  theme_minimal(base_size = 11) +
  theme(strip.placement = "outside",
        # bigger env. variable row labels + bigger x-axis title
        strip.text.y.left = element_text(angle = 0, hjust = 1, size = 11),
        axis.text.x = element_text(size = 11),
        axis.title.x = element_text(size = 13),
        legend.title = element_text(size = 11),
        legend.text = element_text(size = 10),
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.box = "vertical",
        legend.box.just = "center",
        legend.box.spacing = unit(8, "pt"),
        legend.spacing.y = unit(1, "pt"),
        legend.margin = margin(t = 15, b = 0),
        plot.margin = margin(t = 5.5, r = 15, b = 5.5, l = 5.5)) +
  row_grid_theme

#-----------------------------------------------------------------------------
## 2.3) plot panel right: GO enrichment heatmap (BP/CC/MF)
#-----------------------------------------------------------------------------

# create heatmap plot for GO -> columns split by ontology (BP/CC/MF), rows by env_var as before
panel_c_plot <-
  ggplot(panelC_data, aes(x = go_term, y = "row", fill = fdr)) +
  geom_tile(colour = "grey85") +
  facet_grid(rows = vars(env_var), cols = vars(ontology), switch = "both", drop = FALSE,
             scales = "free_x", space = "free_x",
             labeller = labeller(env_var = strip_env_prefix)) +
  scale_fill_gradientn(colours = c("#03045E", "#0077B6", "#00B4D8", "#90E0EF", "#CAF0F8"),
                       limits = c(0, 0.05), name = "GO q-value") +
  # q-values are Benjamini-Hochberg FDR-corrected, from SetRank
  labs(x = "GO term, grouped by ontology (BP = biological process, CC = cellular component, MF = molecular function)",
       y = NULL,
       title = "Gene Ontology Terms across Environmental Variables for PicMin Genes (q-value < 0.25)") +
  theme_minimal(base_size = 11) +
  # own y-axis env. variable, own x-axis strip per ontology (BP/CC/MF)
  theme(strip.placement = "outside",
        strip.text.y.left = element_text(angle = 0, hjust = 1),
        strip.text.x.bottom = element_text(face = "bold"),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 9),
        panel.grid.major.x = element_line(colour = "grey85", linewidth = 0.3),
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank()) +
  row_grid_theme

#-----------------------------------------------------------------------------
## 2.4) combine panels a + b into main summary plot (patchwork, side by side)
#-----------------------------------------------------------------------------
# -> panel c (GO) is now kept out of this combined figure, see 2.5 below

# tag panels a/b directly via labs(tag = ...)
# panel_a_plot <- panel_a_plot + labs(tag = "a")
# panel_b_plot <- panel_b_plot + labs(tag = "b")

# create main plot: panel a + panel b, next to each other
final_plot_main <-
  (panel_a_plot | panel_b_plot) +
  patchwork::plot_layout(widths = c(1, 1)) +
  patchwork::plot_annotation(
    title    = "Cross-Species PicMin Analysis: Signals of Repeated Adaptation",
    theme = theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 16)
    )
  )

# display
final_plot_main
# save: 1200 x 800

#-----------------------------------------------------------------------------
## 2.5) panel c (GO) as standalone supplementary plot
#-----------------------------------------------------------------------------

# -> no tag letter, this is a single-panel supplementary figure on its own

# rename panel c plot
final_plot_supp <- panel_c_plot

# display
final_plot_supp
# save: 800 x 500

#==============================================================================
## 3) summary statistics for thesis text (median / IQR)
#==============================================================================

#-----------------------------------------------------------------------------
## 3.1) median/IQR PicMin gene hit count per env. variable, by method
#-----------------------------------------------------------------------------

# median/IQR for PicMin genes (q-value < 0.25) across env. variables for LFMM and CA
method_summary <- bar_totals %>%
  dplyr::group_by(method) %>%
  dplyr::summarise(
    median_n_genes = round(median(total), 1),
    IQR_n_genes    = round(IQR(total), 1),
    .groups = "drop"
  )

# inspect
method_summary

#-----------------------------------------------------------------------------
## 3.2) save
#-----------------------------------------------------------------------------

# save summary stats and also total bar counts (genes with q-value < 0.25 and counts per env. variable)
write.csv(method_summary, file = "outputs/picmin_4s/picmin_method_summary.csv", row.names = FALSE)
write.csv(bar_totals,     file = "outputs/picmin_4s/picmin_env_method_totals.csv", row.names = FALSE)

