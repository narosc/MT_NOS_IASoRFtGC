#=============================================================================
### 12 - create global map showing all samples
#=============================================================================

# This script creates a global map showcasing all samples for the thesis.
# Input:
# 1) all run tables containing the sample IDs and geo locations
# Output:
# 1. map showing the global distribution of samples

# set correct working directory:
setwd("C:/Users/Admin/Documents/Uni_Zuerich/Masters_Thesis")

# Author: Naroa Olivia Schweizer
# last update: 26.08.2026

# load libraries
library(dplyr)
library(ggplot2)
library(rnaturalearth)
library(ggspatial)
library(ggrepel)
library(patchwork)
library(ggtext)

#-----------------------------------------------------------------------------
## define paths and load files
#-----------------------------------------------------------------------------

# define paths
ptevol_path <- "SraRunTable_PRJNA510810_PteroisVolitans_GeoInfo_coh_addedGeoInf.csv"
amphibic_path <- "SraRunTable_PRJNA294760_Amphiprionbicinctus.csv"
siphtub_path <- "SraRunTable_PRJNA385011_Siphamiatubifer.csv"
epistri_path <- "SraRunTable_PRJEB36904_Ephinephelusstriatus.csv"

# read run tables
ptevol_runtable <- read.csv(ptevol_path)
amphibic_runtable <- read.csv(amphibic_path)
siphtub_runtable <- read.csv(siphtub_path)
epistri_runtable <- read.csv(epistri_path)

## define AOIs spiecis specific
aois <- list(
  "Pterois volitans" = list(
    xlims = c(-90, -60),
    ylims = c(10, 30)
  ),
  "Amphiprion bicinctus" = list(
    xlims = c(30, 60),
    ylims = c(10, 35)
  ),
  "Siphamia tubifer" = list(
    xlims = c(126.75, 128.5),
    ylims = c(26, 27)
  ),
  "Epinephelus striatus" = list(
    xlims = c(-79, -72),
    ylims = c(20.5, 27.5)
  )
)

# equalise plot ratio for all inset maps
equalize_aoi_ratio <- function(aois, target_ratio) {
  lapply(aois, function(a) {
    w <- diff(a$xlims)
    h <- diff(a$ylims)
    cx <- mean(a$xlims)
    cy <- mean(a$ylims)
    current_ratio <- w / h
    
    if (current_ratio < target_ratio) {
      w <- target_ratio * h
    } else if (current_ratio > target_ratio) {
      h <- w / target_ratio
    }
    
    list(
      xlims = c(cx - w / 2, cx + w / 2),
      ylims = c(cy - h / 2, cy + h / 2)
    )
  })
}

# equalise AOIs so ratio in all species the same
aois_eq <- equalize_aoi_ratio(aois, target_ratio = 1.5)

# colour for each species
species_colours <- c(
  "Pterois volitans"     = "#8FC9C2",
  "Amphiprion bicinctus"  = "#FF7F6A",
  "Siphamia tubifer"      = "#6CB6B6",
  "Epinephelus striatus"  = "#F46D75"
)

# legend order
species_legend_order <- c(
  "Epinephelus striatus",
  "Pterois volitans",
  "Amphiprion bicinctus",
  "Siphamia tubifer"
)

#=============================================================================
## 1) prepare input for map
#=============================================================================

#-----------------------------------------------------------------------------
## 1.1) load background map data
#-----------------------------------------------------------------------------

# load base map content
world <- ne_countries(
  scale = "medium",
  returnclass = "sf"
)

# data frame of the species AOI rectangles
min_box_deg <- 6
aoi_rects <- bind_rows(lapply(names(aois), function(sp) {
  xlims <- aois[[sp]]$xlims
  ylims <- aois[[sp]]$ylims
  cx <- mean(xlims); cy <- mean(ylims)
  w <- max(diff(xlims), min_box_deg)
  h <- max(diff(ylims), min_box_deg)
  data.frame(
    species = sp,
    xmin = cx - w / 2,
    xmax = cx + w / 2,
    ymin = cy - h / 2,
    ymax = cy + h / 2
  )
})) |>
  mutate(
    species = factor(species, levels = species_legend_order)
  )

#-----------------------------------------------------------------------------
## 1.2) create function for inset map
#-----------------------------------------------------------------------------

# clean sampling-location names
clean_site_name <- function(x) {
  x <- sub("\\s*\\(\\d{4}\\)\\s*$", "", x)   # drop trailing " (YYYY)"
  x <- sub("^[^:]+:\\s*", "", x)             # drop leading "country: " prefix
  x <- trimws(x)
  gsub("\\s+", " ", x)
}

# function to make species maps
make_inset <- function(data, xlims, ylims, title, point_colour, force_right = character(0)) {
  
  data <- data |>
    mutate(
      geo_loc_name = clean_site_name(geo_loc_name)
    )
  # one dot per distinct sampling coordinate
  sites <- data |>
    distinct(latitude, longitude)
  # sampling site lables -> retrieve one single lables
  site_labels <- data |>
    distinct(geo_loc_name, latitude, longitude) |>
    group_by(geo_loc_name) |>
    mutate(
      lon_c = mean(longitude, na.rm = TRUE),
      lat_c = mean(latitude, na.rm = TRUE),
      dist_c = sqrt((longitude - lon_c)^2 + (latitude - lat_c)^2)
    ) |>
    slice_min(dist_c, n = 1, with_ties = FALSE) |>
    ungroup() |>
    select(geo_loc_name, longitude, latitude)
  site_labels_right <- site_labels |>
    filter(geo_loc_name %in% force_right)
  site_labels_free <- site_labels |>
    filter(!geo_loc_name %in% force_right)
  
  ggplot() +
    geom_sf(
      data = world,
      fill = "grey55",
      colour = "grey40",
      linewidth = 0.25
    ) +
    # sampling locations (species colour)
    geom_point(
      data = sites,
      aes(
        x = longitude,
        y = latitude
      ),
      colour = point_colour,
      size = 2.4,
      alpha = 0.9
    ) +
    # plot one single lable per sampling site lable
    (if (nrow(site_labels_free) > 0) geom_label_repel(
      data = site_labels_free,
      aes(
        x = longitude,
        y = latitude,
        label = geo_loc_name
      ),
      size = 5.2,
      colour = "black",
      fill = alpha("white", 0.85),
      label.size = 0.15,
      label.padding = grid::unit(0.18, "lines"),
      box.padding = 0.55,
      point.padding = 0.25,
      segment.size = 0.3,
      segment.colour = "grey30",
      min.segment.length = 0,
      max.overlaps = Inf,
      force = 2.4,
      force_pull = 0.5,
      seed = 1
    ) else NULL) +
    # labels pinned to the right of their point (force_right)
    (if (nrow(site_labels_right) > 0) geom_label_repel(
      data = site_labels_right,
      aes(
        x = longitude,
        y = latitude,
        label = geo_loc_name
      ),
      size = 5.2,
      colour = "black",
      fill = alpha("white", 0.85),
      label.size = 0.15,
      label.padding = grid::unit(0.18, "lines"),
      box.padding = 0.55,
      point.padding = 0.25,
      segment.size = 0.3,
      segment.colour = "grey30",
      min.segment.length = 0,
      max.overlaps = Inf,
      direction = "y",
      hjust = 0,
      nudge_x = diff(xlims) * 0.14,
      force = 1,
      force_pull = 0.5,
      seed = 1
    ) else NULL) +
    scale_x_continuous(
      breaks = pretty(xlims, n = 4)
    ) +
    scale_y_continuous(
      breaks = pretty(ylims, n = 4)
    ) +
    coord_sf(
      xlim = xlims,
      ylim = ylims,
      expand = FALSE
    ) +
    # scale bar
    annotation_scale(
      location = "bl",
      width_hint = 0.18,
      style = "bar",
      height = grid::unit(0.4, "cm"),
      unit_category = "metric",
      text_cex = 0.9
    ) +
    # north arrow
    annotation_north_arrow(
      location = "tr",
      which_north = "true",
      height = grid::unit(1.1, "cm"),
      width = grid::unit(1.1, "cm"),
      style = north_arrow_fancy_orienteering(text_size = 10)
    ) +
    # species name title
    ggtitle(title) +
    labs(
      x = "Longitude",
      y = "Latitude"
    ) +
    theme_minimal(base_size = 10) +
    theme(
      legend.position = "none",
      # background behind title text -> use ggtext
      plot.title = ggtext::element_textbox(
        face = "bold.italic",
        size = 16,
        colour = "black",
        hjust = 0.5,
        halign = 0.5,
        fill = alpha("white", 0.65),
        box.color = NA,
        padding = margin(4, 10, 4, 10),
        margin = margin(t = 12, b = 10)
      ),
      # axis titles
      axis.title = element_text(size = 14),
      axis.text.x = element_text(
        size = 13,
        colour = "black"
      ),
      axis.text.y = element_text(
        size = 13,
        colour = "black"
      ),
      axis.ticks = element_line(
        colour = "black",
        linewidth = 0.25
      ),
      axis.ticks.length = grid::unit(1.5, "mm"),
      panel.grid.major = element_line(
        colour = "grey85",
        linewidth = 0.2
      ),
      # create colour border
      panel.border = element_blank(),
      panel.background = element_rect(
        colour = NA,
        fill = "white"
      ),
      plot.background = element_rect(
        colour = point_colour,
        fill = "white",
        linewidth = 2.4
      ),
      # create plot margin
      plot.margin = margin(
        t = 16,
        r = 14,
        b = 14,
        l = 14,
        unit = "pt"
      )
    )
}

#=============================================================================
## 2) create map(s)
#=============================================================================

#-----------------------------------------------------------------------------
## 2.1) create inset species specific maps
#-----------------------------------------------------------------------------

# Pterois volitans
p_ptevol <- make_inset(
  data = ptevol_runtable,
  xlims = aois_eq[["Pterois volitans"]]$xlims,
  ylims = aois_eq[["Pterois volitans"]]$ylims,
  title = "Pterois volitans",
  point_colour = species_colours[["Pterois volitans"]]
)

# Amphiprion bicinctus
p_amphibic <- make_inset(
  data = amphibic_runtable,
  xlims = aois_eq[["Amphiprion bicinctus"]]$xlims,
  ylims = aois_eq[["Amphiprion bicinctus"]]$ylims,
  title = "Amphiprion bicinctus",
  point_colour = species_colours[["Amphiprion bicinctus"]],
  force_right = "Djibouti"
)

# Siphamia tubifer
p_siphtub <- make_inset(
  data = siphtub_runtable,
  xlims = aois_eq[["Siphamia tubifer"]]$xlims,
  ylims = aois_eq[["Siphamia tubifer"]]$ylims,
  title = "Siphamia tubifer",
  point_colour = species_colours[["Siphamia tubifer"]]
)

# Epinepehlus striatus
p_epistri <- make_inset(
  data = epistri_runtable,
  xlims = aois_eq[["Epinephelus striatus"]]$xlims,
  ylims = aois_eq[["Epinephelus striatus"]]$ylims,
  title = "Epinephelus striatus",
  point_colour = species_colours[["Epinephelus striatus"]]
)

#-----------------------------------------------------------------------------
## 2.2) create global map
#-----------------------------------------------------------------------------

# create global plot
p_global <- ggplot() +
  geom_sf(
    data = world,
    fill = "grey55",
    colour = "grey40",
    linewidth = 0.2
  ) +
  geom_rect(
    data = aoi_rects,
    aes(
      xmin = xmin,
      xmax = xmax,
      ymin = ymin,
      ymax = ymax,
      colour = species
    ),
    fill = NA,
    linewidth = 1.1
  ) +
  scale_colour_manual(
    values = species_colours,
    breaks = species_legend_order,
    limits = species_legend_order,
    drop = FALSE,
    # legend guide
    guide = guide_legend(
      nrow = 1,
      byrow = TRUE,
      override.aes = list(fill = NA, linewidth = 2.2)
    )
  ) +
  coord_sf(
    xlim = c(-180, 180),
    ylim = c(-60, 80),
    expand = FALSE
  ) +
  labs(
    x = "Longitude",
    y = "Latitude",
    colour = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major = element_line(
      colour = "grey90",
      linewidth = 0.2
    ),
    # axis title
    axis.title = element_text(size = 17),
    axis.text = element_text(size = 12),
    # global legend
    legend.position = c(0.5, 0.07),
    legend.title = element_blank(),
    legend.text = element_text(
      size = 19,
      face = "italic"
    ),
    legend.background = element_rect(
      fill = alpha("white", 0.85),
      colour = "grey40",
      linewidth = 0.3
    ),
    legend.margin = margin(t = 8, r = 14, b = 8, l = 14),
    legend.key.size = grid::unit(22, "pt"),
    legend.spacing.x = grid::unit(6, "mm"),
    plot.margin = margin(t = 5, r = 5, b = 5, l = 5)
  )

#-----------------------------------------------------------------------------
## 2.3) create final figure, combine plots
#-----------------------------------------------------------------------------

# create layout design
layout_design <- "
AB
CC
DE
"

# exact-fit for species maps
true_ratio <- function(xlims, ylims) {
  (diff(xlims) * cos(mean(ylims) * pi / 180)) / diff(ylims)
}

fig_width_in <- 15  # width chosen directly; height is derived below

# each species row has 2 columns -> each inset's cell is half the fig width
species_cell_w_in <- fig_width_in / 2

row1_true_ratio <- mean(c(
  true_ratio(aois_eq[["Epinephelus striatus"]]$xlims, aois_eq[["Epinephelus striatus"]]$ylims),
  true_ratio(aois_eq[["Siphamia tubifer"]]$xlims,     aois_eq[["Siphamia tubifer"]]$ylims)
))
row3_true_ratio <- mean(c(
  true_ratio(aois_eq[["Pterois volitans"]]$xlims,     aois_eq[["Pterois volitans"]]$ylims),
  true_ratio(aois_eq[["Amphiprion bicinctus"]]$xlims, aois_eq[["Amphiprion bicinctus"]]$ylims)
))

row1_cell_h_in <- species_cell_w_in / row1_true_ratio
row3_cell_h_in <- species_cell_w_in / row3_true_ratio

# define global world extent for background map
world_xlim <- c(-180, 180)
world_ylim <- c(-60, 80)
world_ratio <- true_ratio(world_xlim, world_ylim)
global_cell_h_in <- fig_width_in / world_ratio

row_heights <- c(row1_cell_h_in, global_cell_h_in, row3_cell_h_in)
fig_height_in <- sum(row_heights)

# wrap plots
final_map <- wrap_plots(
  A = free(p_epistri),
  B = free(p_siphtub),
  C = free(p_global, side = "tlr"),
  D = free(p_ptevol),
  E = free(p_amphibic),
  design = layout_design
)

# add overall figure title
final_map <- final_map +
  plot_layout(heights = row_heights) +
  plot_annotation(
    title = "Global Distribution of Reef Fish Species and their Sampling Sites",
    theme = theme(
      plot.title = element_text(
        face = "bold",
        size = 24,
        hjust = 0.5,
        margin = margin(b = 10)
      ),
      plot.margin = margin(t = 10, r = 12, b = 6, l = 12)
    )
  )

# save figure
ggsave(
  "global_map_all_samples.png",
  plot = final_map,
  width = fig_width_in,
  height = fig_height_in,
  units = "in",
  dpi = 300,
  bg = "white"
)
