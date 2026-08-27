#=============================================================================
### fourth analysis - compute population structure for samples
#=============================================================================

# This script calculates pairwise Fst (using hierfstat -> Weir and Cockerham 
# (1984) method), where the clusters are computed using the geolocations of 
# each sample with applying a 5 km buffer. Samples within 5 km of each other 
# are grouped into the same cluster.
# Inputs:
# 1) genlight object with imputed/filtered genotype (gt) matrix
# 2) adjusted run table with metadata for samples
# Outputs: 
# 1. a csv file of the pairwise Fst results (because it takes forever to load)
# 2. a map of the samples coloured by regional cluster (input for Fst calculation)
# 3. a heat map of pairwise Fst as a plot

# Author: Naroa Olivia Schweizer
# last update: 04.01.26

# set correct working directory: setwd("C:/Users/Admin/Documents/Uni_Zuerich/Masters_Thesis/CoralReef_fish_local")


# load libraries
library(ggspatial)      # adds nice map annotations (optional)
library(adegenet)       # PCA & multivariate analysis of genetic data
library(dplyr)          # data wrangling (filter, summarize, join)
library(tidyr)          # reshape and tidy data frames
library(hierfstat)      # calculates Fst
library(geosphere)      # calculates buffer for geolocations
library(viridis)        # colour palattes (heatmaps)
library(rnaturalearth)  # for clean world map polygons
library(sf)             # needed for the scale bar and map

#-----------------------------------------------------------------------------
## define paths and load files
#-----------------------------------------------------------------------------

## filtered genligth object (biallelic and sites with less than 10% missing genotypes)
# define path to filtered and imputed genlight object
genlight_obj_path <- "intermed_outputs/epinephelus_striatus/var/estriatus_genlight_filtered.rds"

# read genlight object
genlight_obj <- readRDS(genlight_obj_path)

## Run Table with all infos of samples
# define path to run table
RunTable_path <- "data/EpiStri_RunTable.csv"

# read run table
RunTable <- read.csv(RunTable_path, header = TRUE, stringsAsFactors = FALSE)

## create output directory and define path for image outputs 
# create directory
dir.create('outputs/epinephelus_striatus/pop_struct', showWarnings = FALSE)

# define output path to store visual outputs
output_path <- "outputs/epinephelus_striatus/pop_struct"

## define spatial extent of maps (Area Of Interest: AOI)
# define AOI
xlims = c(-79, -72) # longitude
ylims = c(20.5, 27.5) # latitude

#-----------------------------------------------------------------------------
## 1) prepare gt data
#-----------------------------------------------------------------------------

# extract gt matrix from genlight (0/1/2)  
geno_matrix <- as.matrix(genlight_obj)

# make sure it is only int
geno_matrix <- apply(geno_matrix, 2, as.integer)

#-----------------------------------------------------------------------------
## 2) calculate clusters on geolocation (buffer = 5km), define pop vector
#-----------------------------------------------------------------------------

# make separate table for latitude and longitude to then create distance matrix
coords <- RunTable[, c("longitude", "latitude")]
coords <- as.data.frame(lapply(coords, as.numeric))
rownames(coords) <- RunTable$Run

# compute distance matrix in km
dist_matrix <- distm(coords, fun = distHaversine) / 1000  # distm in m -> divide by 1000 for km

#-------------
# test -> run only within the red sea
#-------------

# # create sf points for coordinates
# coords_sf <- st_as_sf(
#   coords,
#   coords = c("longitude", "latitude"),
#   crs = 4326
# )
# 
# # create a boundry box for the red sea
# redsea_bbox <- st_bbox(
#   c(
#     xmin = 32,  # southern Egypt / Sudan
#     xmax = 44,  # western Arabian coast
#     ymin = 12,  # Eritrea
#     ymax = 30   # Gulf of Suez
#   ),
#   crs = st_crs(4326)
# )
# 
# # make a polygone out of the bbox
# redsea_poly <- st_as_sfc(redsea_bbox)
# 
# # create a mask of the sf points within the readsea polygone
# inside_redsea <- st_within(coords_sf, redsea_poly, sparse = FALSE)[,1]
# 
# # filter out the coordinates that are within the mask
# coords_redsea <- coords[inside_redsea, ]
# 
# # compute distance matrix in km
# dist_matrix <- distm(coords_redsea, fun = distHaversine) / 1000  # distm in m -> divide by 1000 for km
# 
# # add rownames of the samples to the geno_matrix
# rownames(geno_matrix) <- RunTable$Run
# 
# # only retain the samples in the geno_matrix that are also within the red sea
# geno_matrix <- geno_matrix[
#   rownames(geno_matrix) %in% rownames(coords_redsea),
# ]
# 
# # only retain the samples in the run table that are also within the red sea
# RunTable <- RunTable[
#   RunTable$Run %in% rownames(coords_redsea),
# ]

#-------------
# test -> run only outside of the red sea
#-------------

# # create a boundry box for the samples in socrata, oman and dshibuti
# outside_bbox <- st_bbox(
#   c(
#     xmin = 42.0,  # eastern Djibouti / Gulf of Aden
#     xmax = 60.5,  # eastern Oman / Socotra
#     ymin = 10.5,  # southern Djibouti / Arabian Sea
#     ymax = 18   # south of Red Sea (north Gulf of Aden)
#   ),
#   crs = st_crs(4326)
# )
# 
# # make a polygone out of the bbox
# outside_poly <- st_as_sfc(outside_bbox)
# 
# # create a mask of the sf points within the readsea polygone
# outside <- st_within(coords_sf, outside_poly, sparse = FALSE)[,1]
# 
# # filter out the coordinates that are within the mask
# coords_outside <- coords[outside, ]
# 
# # compute distance matrix in km
# dist_matrix <- distm(coords_outside, fun = distHaversine) / 1000  # distm in m -> divide by 1000 for km
# 
# # add rownames of the samples to the geno_matrix
# rownames(geno_matrix) <- RunTable$Run
# 
# # only retain the samples in the geno_matrix that are also within the red sea
# geno_matrix <- geno_matrix[
#   rownames(geno_matrix) %in% rownames(coords_outside),
# ]
# 
# # only retain the samples in the run table that are also within the red sea
# RunTable <- RunTable[
#   RunTable$Run %in% rownames(coords_outside),
# ]

#-------------
# test done -> from here on the regular script
#-------------

# hierarchical clustering
hc <- hclust(as.dist(dist_matrix), method = "average")

# cut tree at 5 km to define clusters
clusters <- cutree(hc, h = 5)

# add cluster info to your metadata
RunTable$Cluster <- clusters

# define population vector from clusters above in geno_all table
pop <- factor(RunTable$Cluster) 

#-----------------------------------------------------------------------------
## 3) prepare gt and pop vector data -> then calculate pairwise Fst
#-----------------------------------------------------------------------------

# final hierfstat data frame (pupulation vector needs to be column 1)
geno_hier <- data.frame(pop = pop, geno_matrix)

# calculate pairwise Fst using weir & cockerham method
fst_pairwise <- pairwise.WCfst(geno_hier)

## save pairwise Fst in table formate!
# covert fst into table format
fst_table <- as.data.frame(fst_pairwise)
# save Fst values as csv (just for saftey as Fst calculation run for ever)
write.csv(
  fst_table,
  file = file.path(output_path, "EpiStri_pairwise_Fst.csv"),
  row.names = FALSE
)

#-----------------------------------------------------------------------------
## 4) visualise data - clusters on map, Fst on map and in table
#-----------------------------------------------------------------------------

## map of the clusters
# get world map polygons
world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")

# create a map coloured by cluster
ggplot() +
  geom_sf(
    data = world,
    fill = "grey95",
    colour = "grey70",
    linewidth = 0.3
  ) +
  geom_point(
    data = RunTable,
    aes(
      x = longitude,
      y = latitude,
      colour = factor(Cluster)
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
    title = "Epinephelus striatus Sampling Sites 
coloured by Regional Cluster [5km]"
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

## visualise a table of the pairwise fst values
# restructure table format to visualse
fst_long <- fst_table %>%
  as.data.frame() %>%
  rownames_to_column(var = "row") %>%
  pivot_longer(
    cols = -row,
    names_to = "col",
    values_to = "fst"
  )

# make sure columns are orderd numerically (1 2 3 4 etc. not 1 10 11 2 3)
fst_long <- fst_long %>%
  mutate(
    row = factor(row, levels = sort(unique(as.numeric(row)))),
    col = factor(col, levels = sort(unique(as.numeric(col))))
  )

# create heatmap of Fst values as image
ggplot(fst_long, aes(x = col, y = row, fill = fst)) +
  geom_tile() +
  scale_fill_viridis_c(option = "C", na.value = "white") +
  coord_fixed() +
  scale_y_discrete(limits = rev) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    plot.title = element_text(hjust = 0.5)
  ) +
  labs(
    title = "Pairwise Fst of Regional Clusters for Epinephelus striatus",
    x = "Regional Cluster",
    y = "Regional Cluster",
    fill = "Fst"
  )
