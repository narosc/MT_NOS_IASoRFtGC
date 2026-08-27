#=============================================================================
### fifth analysis - analyse environmental data - animations, PCA
#=============================================================================

# This script analyses the environmental data, creates animated heat maps of 
# different time series variables and calculates PCA on the environmental 
# variables.
# Inputs:
# 1) meta data for the RECIFs data (not needed for calculation, just understanding)
# 2) RECIFs_anim -> monthly averages per year for animations trimmed to convex hull
# 3) RECIFs_samples -> avg & sd for sample reef cells
# 4) adjusted run table with metadata for samples
# Outputs: 
# 1. a csv file of environmental PCA scores with coordinates & sample names
# 2. an R object of the environmental PCA results (total, takes long to load)
# 3. a scree plot with broken stick criterion of PCA axis
# 4. two bi plots (avg & sd) coloured by sampling region 
# 5. a correlation plot -> pairwise Pearson's R for all environmental variables
# (6. a map of sampling sites coloured by sampling region) -> already created

# It creates two outputs: 
# 1) Animations of the longer time series as a visual output to look at the 
# environmental change over the years.
# 2) A patchwork plot of PCA results as well as the saved data frame and PCA 
# result (R Object) for the PCA for the further analysis (LFMM).

# Note: Check if CRS in both the sample coordinates and the environmental data
# are the same.

# Author: Naroa Olivia Schweizer
# last update: 25.01.2026

# set correct working directory: setwd("C:/Users/Admin/Documents/Uni_Zuerich/Masters_Thesis/CoralReef_fish_local")

# load libraries
library(ggspatial)      # adds nice map annotations (optional)
library(ggplot2)        # general-purpose plotting
library(dplyr)          # data wrangling (filter, summarize, join)
library(tidyr)          # reshape and tidy data frames
library(viridis)        # colour palattes (heatmaps)
library(grid)           # needed to compute scale bar
library(magick)         # create animations (for time series)
library(ggnewscale)     # for multiple color scales in one plot
library(ggcorrplot)     # plot correlation graphs
library(rnaturalearth)  # for clean world map polygons
library(ggdendro)       # for creating dendro plots

#-----------------------------------------------------------------------------
## define paths and load files
#-----------------------------------------------------------------------------

## sample
# define path to run table
RunTable_path <-"data/EpiStri_RunTable.csv"

# read run table
RunTable <- read.csv(RunTable_path, header = TRUE, stringsAsFactors = FALSE)

## RECIFs
# define path to meta data for recifs
meta_recifs_path <- "intermed_outputs/RECIFs_inputs/meta_variables_recifs.csv"

# read meta data recifs
meta_recifs <- read.csv(meta_recifs_path)

# define path to filtered recifs data for animations
RECIFs_anim_path <- "intermed_outputs/epinephelus_striatus/env_data_inputs/EpiStri_RECIFs_anim.csv"

# read filtered recifs data for animations
RECIFs_anim <- read.csv(RECIFs_anim_path)

# define path to filtered recifs data for PCA (avg and sd)
RECIFs_sample_path <- "intermed_outputs/epinephelus_striatus/env_data_inputs/EpiStri_RECIFs_sample.csv"

# read filtered recifs data for PCA (avg and sd data for begining to 2011)
RECIFs_sample <- read.csv(RECIFs_sample_path, header = TRUE)
# add sample row names to env. data sample
rownames(RECIFs_sample) <- RunTable$Run

## outputs
# create directory
dir.create('outputs/epinephelus_striatus/animations', showWarnings = FALSE)
# define output path to store animation outputs
output_path <- "outputs/epinephelus_striatus/animations"

# define output path to save PCA outputs
output_path_PCA <- "outputs/epinephelus_striatus"

# define output path to save preped env data
save_path <- "intermed_outputs/epinephelus_striatus/env_data_inputs/"

## define spatial extent of maps (Area Of Interest: AOI)
# define AOI
xlims = c(-79, -72) # longitude
ylims = c(20.5, 27.5) # latitude

## define species specific end year for animation calculations
# define Yend
yend <- 2014

#=============================================================================
## 1) visualise heat maps of the data - create animations 
#=============================================================================

#-----------------------------------------------------------------------------
## 1.1) visualise heat maps of the data 
#-----------------------------------------------------------------------------

## load base map, define unchangeable inputs
# get world map polygons
world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")

## define changeable variables (choose colour palatte for each variable)
variable_info <- list(
  CHL = list(years=1997:yend, name="Avg CHL [mg m-3]", palette="viridis"),
  DHW = list(years=1985:yend, name="Avg Degree Heating Week [°C-week]", palette="plasma"),
  FE  = list(years=1993:yend, name="Avg Iron [mmol m-3]", palette="cividis"),
  O2  = list(years=1993:yend, name="Avg O2 [mmol m-3]", palette="cividis"),
  PH  = list(years=1993:yend, name="Avg pH", palette="viridis"),
  NO3 = list(years=1993:yend, name="Avg Nitrate [mmol m-3]", palette="viridis"),
  PO4 = list(years=1993:yend, name="Avg Phosphate [mmol m-3]", palette="viridis"),
  SCV = list(years=1993:yend, name="Avg Sea Water Velocity [m/s]", palette="viridis"),
  SPM = list(years=1997:yend, name="Avg Suspended Matter [g/m3]", palette="viridis"),
  SSS = list(years=1993:yend, name="Avg Sea Surface Salinity [1e-3]", palette="viridis"),
  SST = list(years=1985:yend, name="Avg Sea Surface Temperature [°C]", palette="magma")
)

#-----------------------------------------------------------------------------
## 1.2) create the base map for the animations
#-----------------------------------------------------------------------------

# base map shared by all frames
base_map <- ggplot() +
  geom_sf(
    data = world,
    fill = "grey95",
    colour = "grey70",
    linewidth = 0.3
  ) +
  coord_sf(
    xlim = xlims,
    ylim = ylims,
    expand = FALSE
  ) +
  annotation_scale(
    location      = "bl",
    width_hint    = 0.2,
    style         = "bar",
    bar_cols      = c("black", "white"),
    height        = unit(0.4, "cm"),
    unit_category = "metric",
    dist_unit     = "km"
  ) +
  annotation_north_arrow(
    location    = "tr",
    which_north = "true",
    style       = north_arrow_fancy_orienteering()
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major = element_line(colour = "grey90"),
    legend.position  = "right"
  ) +
  labs(
    x = "Longitude",
    y = "Latitude"
  )

#-----------------------------------------------------------------------------
## 1.3) run loop over all time series variables 
#-----------------------------------------------------------------------------

# create and save animations for a given variable using this function
animate_variable <- function(var_code, years, var_name, color_option, output_file) {
  
  var_cols <- paste0("s_allDB_", var_code, "_010.", var_code, "_", years)
  if (!all(var_cols %in% colnames(RECIFs_anim))) {
    stop(paste("Columns missing for", var_code))
  }
  
  var_min <- min(RECIFs_anim[, var_cols], na.rm = TRUE)
  var_max <- max(RECIFs_anim[, var_cols], na.rm = TRUE)
  if (!is.finite(var_min)) var_min <- 0
  if (!is.finite(var_max)) var_max <- 1
  
  png_files <- character(0)
  
  for (i in seq_along(years)) {
    col <- var_cols[i]
    yr  <- years[i]
    
    p <- base_map +
      geom_point(
        data = RECIFs_anim,
        aes(
          x     = allDB_coords.lon,
          y     = allDB_coords.lat,
          color = .data[[col]]
        ),
        size  = 2,
        alpha = 0.8
      ) +
      scale_color_viridis_c(
        option = color_option,
        name   = var_name,
        limits = c(var_min, var_max)
      ) +
      labs(
        title = paste(var_name, yr)
      )
    
    tmp_file <- tempfile(
      pattern = paste0(var_code, "_", yr, "_"),
      fileext = ".png"
    )
    
    ggsave(
      filename = tmp_file,
      plot     = p,
      width    = 9,
      height   = 5,
      dpi      = 150
    )
    
    png_files <- c(png_files, tmp_file)
  }
  
  anim <- image_animate(
    image_join(image_read(png_files)),
    fps = 2
  )
  
  image_write(anim, path = file.path(output_path, output_file))
  unlink(png_files)
  
  invisible(NULL)
}

# loop all individual inputs through the function
for (var in names(variable_info)) {
  info <- variable_info[[var]]
  animate_variable(
    var,
    info$years,
    info$name,
    info$palette,
    paste0(var, "_anim.gif")
  )
}

#-----------------------------------------------------------------------------
## 1.4) general functioning code for spatial heat map
#-----------------------------------------------------------------------------

# plot
ggplot() +
  geom_sf(
    data = world,
    fill = "grey95",
    colour = "grey70",
    linewidth = 0.3
  ) +
  geom_point(
    data = RECIFs_anim,
    aes(
      x = allDB_coords.lon,
      y = allDB_coords.lat,
      color = BATHY
    ),
    size  = 2,
    alpha = 0.8
  ) +
  coord_sf(
    xlim = xlims,
    ylim = ylims,
    expand = FALSE
  ) +
  scale_color_viridis_c(
    option = "viridis",
    name   = "Avg Depth [m]"
  ) +
  annotation_scale(
    location      = "bl",
    width_hint    = 0.2,
    style         = "bar",
    bar_cols      = c("black", "white"),
    height        = unit(0.4, "cm"),
    unit_category = "metric",
    dist_unit     = "km"
  ) +
  annotation_north_arrow(
    location    = "tr",
    which_north = "true",
    style       = north_arrow_fancy_orienteering()
  ) +
  labs(
    x     = "Longitude",
    y     = "Latitude",
    title = "Average Depth [m] for Reef Sites"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major = element_line(colour = "grey90"),
    legend.position  = "right"
  )


#=============================================================================
## 2) run PCA on all variables
#=============================================================================

#-----------------------------------------------------------------------------
## 2.1) run PCA on the geographically filtered data
#-----------------------------------------------------------------------------

# define columns to be excluded (don't want lat and lon to be part of PCA)
exclude_cols <- c("allDB_coords.lon", "allDB_coords.lat", 
                  "sample_lat", "sample_lon")

# select df to calculate PCA on (exclude lat/lon)
pca_input <- RECIFs_sample %>%
  select(
    -all_of(exclude_cols)
  )

# calculate PCA -> scale TRUE (normalize the data each input 0 to 1)
pca_result <- prcomp(pca_input, scale. = TRUE)

# extract PCA scores -> convert to df
pca_scores <- as.data.frame(pca_result$x)

# attach excluded columns to df
pca_scores <- cbind(
  RECIFs_sample[, exclude_cols],
  pca_scores
)

# variance explained by each PC
var_explained <- 100 * (pca_result$sdev^2) / sum(pca_result$sdev^2)

#-----------------------------------------------------------------------------
## 2.2) visualise PCA outputs
#-----------------------------------------------------------------------------

# join run table and pca scores (also add rowname to column)
pca_scores_meta <- pca_scores %>%
  rownames_to_column(var = "sample") %>%
  inner_join(
    RunTable,
    by = c("sample" = "Run")
  )

## calculate broken stick criterion
# number of PCs
K <- length(var_explained)

# broken stick expectations (in proportion)
broken_stick <- sapply(1:K, function(i) {
  sum(1 / (i:K)) / K
})

# convert to percent
broken_stick <- 100 * broken_stick

# create df for scree plot of var explained
scree_df <- data.frame(
  PC = seq_len(K),
  Variance = var_explained,
  BrokenStick = broken_stick
)

# reshape data to long formate
scree_long <- scree_df %>%
  pivot_longer(
    cols = c(Variance, BrokenStick),
    names_to = "Type",
    values_to = "Value"
  )

# plot scree plot
ggplot(scree_long, aes(x = PC, y = Value, colour = Type, linetype = Type)) +
  geom_line(linewidth = 0.8) +
  geom_point(
    data = subset(scree_long, Type == "Variance"),
    size = 2
  ) +
  scale_colour_manual(
    values = c("Variance" = "black", "BrokenStick" = "firebrick"),
    labels = c("Broken stick", "Observed variance")
  ) +
  scale_linetype_manual(
    values = c("Variance" = "solid", "BrokenStick" = "dashed"),
    labels = c("Broken stick", "Observed variance")
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.major = element_line(colour = "grey90"),
    legend.position = c(0.85, 0.85),   # top-right inside plot
    legend.title = element_blank()
  ) +
  labs(
    x = "# PC",
    y = "Variance Explained [%]",
    title = "Scree Plot of Environmental PCA"
  )

## compute biplot for PCA outputs
# extract loadings from PCA result
loadings <- as.data.frame(pca_result$rotation)

# add variable names as column
loadings$varname <- rownames(loadings)

# identify type -> avg or sd
loadings$type <- ifelse(grepl("_sd$", loadings$varname), "sd", "avg")

# extract variable code (CHL, DHW, PO4 etc.)
loadings$varcode <- sub("^s_allDB_([A-Z0-9]+)_(avg|sd)$", "\\1", loadings$varname)

# scale fator -> adjust scale to fit the plot
scale_factor <- 15 
# run the scale factor on the loading matrix
pc_cols <- grep("^PC[0-9]+$", colnames(loadings), value = TRUE) # only scale the PC columns
loadings_scaled <- loadings
loadings_scaled[, pc_cols] <- loadings[, pc_cols] * scale_factor

# split into two data sets 
load_sd  <- loadings_scaled[loadings_scaled$type == "sd", ]
load_avg <- loadings_scaled[loadings_scaled$type == "avg", ]

# offset for labels in the plot (where the VarCode on the arrow)
offset <- 0.4

# biplot for env PCA -> on avg -> datapoints coloured by sampling region
ggplot() +
  geom_point(
    data = pca_scores_meta,
    aes(x = PC1, y = PC2, color = geo_loc_name),
    size = 2,
    alpha = 0.8
  ) +
  scale_colour_hue(name = "Sampling Region") +
  ggnewscale::new_scale_color() +
  geom_segment(
    data = load_avg,
    aes(x = 0, y = 0, xend = PC1, yend = PC2),
    color = "grey20",
    arrow = arrow(length = unit(0.25, "cm")),
    linewidth = 0.8
  ) +
  geom_text(
    data = load_avg,
    aes(
      x = PC1 + sign(PC1) * offset,
      y = PC2 + sign(PC2) * offset,
      label = varcode
    ),
    color = "grey20",
    show.legend = FALSE,
    size = 3.5
  ) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Biplot of Environmental Variables (Averages): PC1 vs PC2",
    x = paste0("PC1 (", round(var_explained[1], 1), "%)"),
    y = paste0("PC2 (", round(var_explained[2], 1), "%)")
  ) +
  theme(legend.position = "right")

# biplot for env PCA -> on sd -> datapoints coloured by sampling region
ggplot() +
  geom_point(
    data = pca_scores_meta,
    aes(x = PC1, y = PC2, color = geo_loc_name),
    size = 2,
    alpha = 0.8
  ) +
  scale_colour_hue(name = "Sampling Region") +
  ggnewscale::new_scale_color() +
  geom_segment(
    data = load_sd,
    aes(x = 0, y = 0, xend = PC1, yend = PC2),
    color = "grey20",
    arrow = arrow(length = unit(0.25, "cm")),
    linewidth = 0.8
  ) +
  geom_text(
    data = load_sd,
    aes(
      x = PC1 + sign(PC1) * offset,
      y = PC2 + sign(PC2) * offset,
      label = varcode
    ),
    color = "grey20",
    show.legend = FALSE,
    size = 3.5
  ) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Biplot of Environmental Variables (st. dev.): PC1 vs PC2",
    x = paste0("PC1 (", round(var_explained[1], 1), "%)"),
    y = paste0("PC2 (", round(var_explained[2], 1), "%)")
  ) +
  theme(legend.position = "right")

## plot a map of sampling sites (coloured by region)
# get world map polygons -> already loaded above
# create a map for the sampling sites
ggplot() +
  geom_sf(
    data = world,
    fill = "grey95",
    colour = "grey70",
    linewidth = 0.3
  ) +
  geom_point(
    data = pca_scores_meta,
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

#-----------------------------------------------------------------------------
## 2.3) save PCA score matrix, prept env variables for species of interest
#-----------------------------------------------------------------------------

# save PCA scores as csv
save_file <- file.path(save_path, "env_PCA_scores.csv")
write.csv(pca_scores, save_file, row.names = TRUE)

# save PCA results (list of different objects) an r object
save(pca_result,
     file = file.path(save_path, "env_PCA_results.RData")
)

#=============================================================================
## 3) Check env. variables for correlation/multicollinearity  
#=============================================================================

# compute pairwise correlations using pearson's R
env_cor <- cor(pca_input, use = "pairwise.complete.obs", method = "pearson")

# plot correlation plot
ggcorrplot(
  env_cor,
  type = "lower",
  lab = TRUE,        # print Pearson's r
  lab_size = 1,     # size of r values inside squares
) +
  scale_fill_gradient2(name = "Pearson's R") +
  ggtitle("Correlation matrix of env var for Epinephelus striatus") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 12),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 6),
    axis.text.y = element_text(size = 6)
  )

#-----------------------------------------------------------------------------
## 3.1) calculate correlation groups of env. variables
#-----------------------------------------------------------------------------

# convert correlation to distance
cor_dist <- as.dist(1 - abs(env_cor))

# hierarchical clustering
hc_env <- hclust(cor_dist, method = "average") 

# plot dedrogram
ggdendro::ggdendrogram(
  hc_env,
  rotate = TRUE,
  size = 0.5
) +
  ggtitle("Hierarchical clustering of environmental variables") +
  ylab("1 − |Pearson's r|") +
  xlab("Env Variable") +
  theme_minimal(base_size = 12)

# group into 5 clusters
var_groups <- cutree(hc_env, k = 5)

# crate df out of var_groups
var_groups_df <- as.data.frame(var_groups)

# prepare data for coloured dendro plot
dend <- ggdendro::dendro_data(hc_env, type = "rectangle")
# extract lables and var_groups from var groups df
label_df <- dend$labels
label_df$group <- factor(var_groups_df[label_df$label, 1])

# identify leaf segments
leaf_segments <- dend$segments[dend$segments$yend == 0, ]
# group leaf segments
leaf_segments$group <- factor(
  label_df$group[match(leaf_segments$xend, label_df$x)]
)

# define a fixed y position for env variables (so that lables are always shown!)
label_df$y_fixed <- -0.01

# choose 5 colors from dark to light red
#my_colors_base <- c("#A33F39", "#C25959", "#9EBFC4", "#68B6BD", "#6CA9B8")
my_colors <- c("#4D96A0", "#86A7AC", "#38C3C9", "#C25959","#A33F39")

# plot dengrogram again coloured by group
ggplot() +
  geom_segment(
    data = dend$segments,
    aes(x = x, y = y, xend = xend, yend = yend),
    colour = "grey70",
    linewidth = 0.5
  ) +
  geom_segment(
    data = leaf_segments,
    aes(x = x, y = y, xend = xend, yend = yend, colour = group),
    linewidth = 1,
    show.legend = FALSE
  ) +
  geom_text(
    data = label_df,
    aes(x = x, y = y_fixed, label = label, colour = group),
    size = 2.5,
    hjust = 1,
    show.legend = FALSE
  ) +
  scale_colour_manual(values = my_colors) +
  coord_flip(clip = "off") +
  labs(
    title = "Hierarchical clustering of environmental variables",
    x = NULL,
    y = "1 − |Pearson's r|"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid   = element_blank(),
    axis.text.y  = element_blank(),
    axis.ticks.y = element_blank(),
    plot.title   = element_text(hjust = 0.5),
    plot.margin  = margin(5, 5, 5, 50)
  )

