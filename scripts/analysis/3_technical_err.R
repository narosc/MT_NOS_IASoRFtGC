#=============================================================================
### thrid analysis - check if technical artefacts correlate with multiqc data
#=============================================================================

# This script looks if the multiqc data correlates with the genotypes by 
# computing a PCoA (Principal Coordinates Analysis) from the filtered genlight 
# object and then calculates the Pearson's R for selected variables.
# Inputs:
# 1) genlight object filtered for relevant samples
# 2) multiqc_data.json file from last multiqc (variant calling) run in initial pipeline
# 3) adjusted run table with metadata for samples
# Outputs: 
# 1. a patchwork plot summarising variance in genotype matrix and correlation 
# to multiqc metrics
# 2. a map showing the sample locations coloured by sampling region (save manually)

# Author: Naroa Olivia Schweizer
# last update: 19.01.26

# set correct working directory: setwd("C:/Users/Admin/Documents/Uni_Zuerich/Masters_Thesis/CoralReef_fish_local")

# load libraries
library(ggspatial)      # adds nice map annotations
library(adegenet)       # PCA & multivariate analysis of genetic data
library(ggplot2)        # general-purpose plotting
library(dplyr)          # data wrangling (filter, summarize, join)
library(TidyMultiqc)    # read and process MultiQC output tables
library(patchwork)      # to frame/combine multiple plots inside one image
library(tibble)         # handling tidyverse-friendly version of df
library(rnaturalearth)  # for clean world map polygons
library(purrr)          # iterating over lists/vectors without writing for loops
library(jsonlite)       # read JSON file for multiqc -> only needed for E. striatus

#-----------------------------------------------------------------------------
## define paths/species specific information and load files
#-----------------------------------------------------------------------------

## filtered genligth object (biallelic and sites with less than 10% missing genotypes)
# define path to filtered and imputed genlight object
genlight_obj_path <- "intermed_outputs/epinephelus_striatus/var/estriatus_genlight_filtered.rds"

# read genlight object
genlight_obj <- readRDS(genlight_obj_path)

## MultiQC file
# define path to multiqc
mqc_path <- "intermed_outputs/epinephelus_striatus/multiqc_LDvc/multiqc_data/multiqc_data.json"

# read MultiQC
mqc <- TidyMultiqc::load_multiqc(mqc_path)
# for epinephelus striatus
mqc_path <- "intermed_outputs/epinephelus_striatus/multiqc_LDvc/multiqc_data/multiqc_data.json"
run_table <- read.csv("data/EpiStri_RunTable.csv")
keep_ids <- run_table$Run

raw <- jsonlite::read_json(mqc_path, simplifyVector = FALSE)

# only keep samples that are actually in your cohort (".coverage" rows share the base ID)
keep_sample <- function(sample_id) sub("\\.coverage$", "", sample_id) %in% keep_ids

is_bad <- function(val) {
  is.null(val) || length(val) != 1 || (is.numeric(val) && is.nan(val))
}

raw$report_general_stats_data <- map(raw$report_general_stats_data, function(module) {
  module <- module[keep(names(module), keep_sample)]
  map(module, function(sample_data) sample_data[!map_lgl(sample_data, is_bad)])
})

tmp_path <- tempfile(fileext = ".json")
jsonlite::write_json(raw, tmp_path, auto_unbox = TRUE)

mqc <- TidyMultiqc::load_multiqc(tmp_path)

## Run Table with all infos of samples
# define path to run table
RunTable_path <- "data/EpiStri_RunTable.csv"

# read run table
RunTable <- read.csv(RunTable_path, header = TRUE, stringsAsFactors = FALSE)

## create output directory and define path for image outputs 
# create directory
dir.create('outputs/epinephelus_striatus', showWarnings = FALSE)
# define output path
output_path <- "outputs/epinephelus_striatus"

## define spatial extent of maps (Area Of Interest: AOI)
# define AOI
xlims = c(-79, -72) # longitude
ylims = c(20.5, 27.5) # latitude

#-----------------------------------------------------------------------------
## 1) extract genotype matrix and calculate PCoA
#-----------------------------------------------------------------------------

# load genlight object as df
gt <- as.data.frame(as.matrix(genlight_obj))

# retrieve distance from genlight object
geno_dist <- dist(gt, method = "euclidean")

# run PCoA on distance matrix of genlight object
pcoa_geno <- cmdscale(
  geno_dist,
  k = 20,
  eig = TRUE
)

# reframe PCoA points into usable df
# -> clean sample names so that file is ready for inner join
geno_table <- as.data.frame(pcoa_geno$points)
# rename axes PCoA
colnames(geno_table) <- paste0("PCo", seq_len(ncol(geno_table)))

#-----------------------------------------------------------------------------
## 2) join all tables - clean table for analysis
#-----------------------------------------------------------------------------

## join with QC metrics
# first remove cohort and cohort_biallelic_fmiss from multiqc table
mqc <- mqc %>%
  filter(!metadata.sample_id %in% "cohort_biallelic_fmiss_maf")
# create subset with only metrics unrelated to coverage
mqc_main <- mqc %>%
  filter(!grepl("\\.coverage$", metadata.sample_id))
# create subset with only metrics related to coverage, mutate metadata.sample_id (delete .coverage suffix)
mqc_cov <- mqc %>%
  filter(grepl("\\.coverage$", metadata.sample_id)) %>%
  mutate(
    metadata.sample_id = sub("\\.coverage$", "", metadata.sample_id)
  )
# select only coverage related columns
mqc_cov <- mqc_cov %>%
  select(
    metadata.sample_id,
    general.numreads,
    general.covbases,
    general.coverage,
    general.meandepth,
    general.meanbaseq,
    general.meanmapq
  )
# merge the two mqc subsets together
mqc_fixed <- mqc_main %>%
  left_join(
    mqc_cov,
    by = "metadata.sample_id",
    suffix = c(".main", ".cov")
  ) %>%
  mutate(
    general.coverage  = coalesce(general.coverage.main,  general.coverage.cov),
    general.covbases  = coalesce(general.covbases.main,  general.covbases.cov),
    general.meandepth = coalesce(general.meandepth.main, general.meandepth.cov)
  ) %>%
  select(
    -ends_with(".main"),
    -ends_with(".cov")
  )
# join PCoA pints of gt matrix to multiqc table (also add rowname to column)
geno_qc <- geno_table %>%
  rownames_to_column(var = "Sample") %>%
  inner_join(
    mqc_fixed,
    by = c("Sample" = "metadata.sample_id")
  )

## join with RunTable (geography, sample location etc.)
geno_all <- geno_qc %>%
  inner_join(RunTable, by = c("Sample" = "Run"))

#=============================================================================
## 3) create output image of different maps/plots
#=============================================================================

#-----------------------------------------------------------------------------
## 3.1) create scree plot and PCoA plots (difference of genetic distance)
#-----------------------------------------------------------------------------

## calculate explained variance of all PCoA axis
var_explained <- 100 * pcoa_geno$eig / sum(pcoa_geno$eig)

# scree plot of PCo #
p1 <- ggplot(data.frame(PCo = seq_along(var_explained), Var = var_explained),
             aes(PCo, Var)) +
  geom_line(linewidth = 0.6, colour = "grey40") +
  geom_point(size = 2, colour = "black") +
  theme_minimal(base_size = 14) +
  theme(panel.grid.major = element_line(colour = "grey90")) +
  labs(
    x = "PCo #",
    y = "Variance Explained [%]",
    title = "Scree Plot for genotype PCoA"
  )

## create genotype PCoA (1,2 | 3,4 | 5,6) plot where colour = sampling region
# define the pairs of PCs you want to plot
pco_pairs <- list(c(1, 2), c(3, 4), c(5, 6))

# create list to store seperate images in...
p_list <- list()

# build plots pc 1,2 | 3,4 | 5,6
for (i in seq_along(pco_pairs)) {
  
  pcos <- pco_pairs[[i]]
  
  p_list[[i]] <- ggplot(
    geno_all,
    aes(
      x = .data[[paste0("PCo", pcos[1])]],
      y = .data[[paste0("PCo", pcos[2])]],
      colour = geo_loc_name
    )
  ) +
    geom_point(size = 2, alpha = 0.8) +
    theme_minimal(base_size = 14) +
    labs(
      x = paste0("PCo", pcos[1], " [", round(var_explained[pcos[1]], 1), "%]"),
      y = paste0("PCo", pcos[2], " [", round(var_explained[pcos[2]], 1), "%]"),
      colour = "Sampling Region"
    )
}

p2 <- p_list[[1]]
p3 <- p_list[[2]]
p4 <- p_list[[3]]

#-----------------------------------------------------------------------------
## 3.2) create plots for multiqc data - genotype correlation
#-----------------------------------------------------------------------------

## create PCo1/2 vs. multiqc metrics plots (calculate pearsons R), sampling region is colour

# define PCos and metrics
pcos     <- c("PCo1", "PCo1", "PCo1", "PCo1", "PCo2", "PCo2", "PCo2", "PCo2")
metrics <- c("general.reads_mapped", "general.percent_gc", "general.total_length", "general.covbases",
             "general.reads_mapped", "general.percent_gc", "general.total_length", "general.covbases")
x_labels <- c("Mapped Reads [bp]", "GC Content [%]", "Read Length [bp]", "Coverage [bp]",
              "Mapped Reads [bp]", "GC Content [%]", "Read Length [bp]", "Coverage [bp]")

# create list to store seperate plots in
qc_plots <- list()

# build plots PCo1/PCo2 vs. MultiQC metrics -> creates p5–p12 automatically
for (i in seq_along(metrics)) {
  
  pco    <- pcos[i]
  metric <- metrics[i]
  
  corV <- cor(
    geno_all[[pco]],
    geno_all[[metric]],
    use = "complete.obs"
  )
  
  qc_plots[[i]] <- ggplot(
    geno_all,
    aes(
      x = .data[[metric]],
      y = .data[[pco]],
      colour = geo_loc_name
    )
  ) +
    geom_point(size = 2, alpha = 0.8) +
    theme_minimal(base_size = 14) +
    theme(
      plot.caption = element_text(
        hjust = 0.5,   # center horizontally
        size = 11,
        face = "bold"
      )
    ) +
    labs(
      x = x_labels[i],
      y = pco,
      colour = "Sampling Region",
      caption = paste0("R = ", round(corV, 2))
    )
}

p5 <- qc_plots[[1]]; p6 <- qc_plots[[2]]
p7 <- qc_plots[[3]]; p8 <- qc_plots[[4]]
p9 <- qc_plots[[5]]; p10 <- qc_plots[[6]]
p11 <- qc_plots[[7]]; p12 <- qc_plots[[8]]

#-----------------------------------------------------------------------------
## 3.3) asseble subplots into one big image
#-----------------------------------------------------------------------------

# row 1: scree plot + PCA plots
row1 <- p1 + p2 + p3 + p4 + 
  plot_layout(guides = "collect", widths = rep(1, 4)) & # equal widths 
  theme(legend.position = "right") & 
  theme(plot.title = element_blank()) # remove individual titles 

# row 2: PC1 vs metrics 
row2 <- p5 + p6 + p7 + p8 + 
  plot_layout(guides = "collect", widths = rep(1, 4)) & 
  theme(legend.position = "right") & 
  theme(plot.title = element_blank()) 

# row 3: PC2 vs metrics 
row3 <- p9 + p10 + p11 + p12 + 
  plot_layout(guides = "collect", widths = rep(1, 4)) & 
  theme(legend.position = "right") & 
  theme(plot.title = element_blank())

## combine rows
final_plot_without <- row1 / row2 / row3 +
  plot_layout(heights = rep(1, 3)) +  # all rows same height
  plot_annotation(
    title = "Technical Error Check - Epinephelus striatus",
    theme = theme(
      plot.title = element_text(
        size = 16,
        face = "bold",
        hjust = 0.5   # centered
      ))
    )

# display
final_plot_without

## save the plot
ggsave(
  filename = file.path(output_path, "epistri_technical_err_PCoA_combined.png"), 
  plot = final_plot_without, 
  width = 16,
  height = 11,  
  dpi = 300
)

#-----------------------------------------------------------------------------
## 3.4) create map of samples
#-----------------------------------------------------------------------------

# get world map polygons
world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")

# create a map for the sampling sites
ggplot() +
  geom_sf(
    data = world,
    fill = "grey95",
    colour = "grey70",
    linewidth = 0.3
  ) +
  geom_point(
    data = geno_all,
    aes(
      x = longitude,
      y = latitude,
      colour = geo_loc_name
    ),
    size = 2,
    alpha = 0.8
  ) +
  coord_sf(
    xlim = xlims,
    ylim = ylims,
    expand = FALSE
  ) +
  scale_colour_hue(name = "Sampling Region") +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major = element_line(colour = "grey90"),
    legend.position = "right"
  ) +
  labs(
    x = "Longitude",
    y = "Latitude",
    title = "Epinephelus striatus Sampling Sites"
  ) +
annotation_scale(
  location     = "bl",         # bottom-left
  width_hint   = 0.2,          # relative width of scale bar
  style    = "bar",
  bar_cols = c("black", "white"),
  height       = unit(0.4, "cm"),
  unit_category = "metric",
  dist_unit    = "km"
) +
annotation_north_arrow(
  location = "tr",
  which_north = "true",
  style = north_arrow_fancy_orienteering()
)


# old version of ggspatial
ggplot() +
  geom_sf(
    data = world,
    fill = "grey95",
    colour = "grey70",
    linewidth = 0.3
  ) +
  geom_point(
    data = geno_all,
    aes(
      x = longitude,
      y = latitude,
      colour = geo_loc_name
    ),
    size = 2,
    alpha = 0.8
  ) +
  coord_sf(
    xlim = xlims,
    ylim = ylims,
    expand = FALSE
  ) +
  scale_colour_hue(name = "Sampling Region") +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major = element_line(colour = "grey90"),
    legend.position = "right"
  ) +
  labs(
    x = "Longitude",
    y = "Latitude",
    title = "Pterois volitans Sampling Sites"
  ) +
  annotation_scalebar(
    location     = "bl",         # bottom-left
    width_hint   = 0.2,          # relative width of scale bar
    style    = "bar",
    bar_cols = c("black", "white"),
    height       = unit(0.4, "cm"),
    unit_category = "metric",
    dist_unit    = "km"
  ) +
  annotation_north_arrow(
    location = "tr",
    which_north = "true",
    style = north_arrow_fancy_orienteering()
  )
