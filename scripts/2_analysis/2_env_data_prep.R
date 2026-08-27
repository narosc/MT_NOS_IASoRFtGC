#=============================================================================
### second analysis - prepare environmental data for analysis
#=============================================================================

# This script creates different subsets of environmental input data from RECIFs.
# Inputs:
# 1) meta data for the RECIFs data (not needed for calculation, just understanding)
# 2) compressed zip file of RECIFs data (-> uncompressed an R-Object)
# 3) adjusted run table with metadata for samples
# General Output: -> this part once ran does not need to be run again! 
# 1. save an R-Object with 10km reef cells worldwide for all variables, where 
# specifically the monthly averages for time series variables is calculated
# Specific Outputs: specific to each species -> needs to be run specifically
# 2. save an R-Object of 10km reef cells worldwide for all variables, where
# the yearly averages and standard deviations within species specific time range
# have been calculated
# 3. saves a csv with all variables (from output 2: avg and sd) filtered to
# nearest neighbour reef cell to sample coordinates 
# 4. saves a csv with all variables (from output 3: monthly averages per year) 
# filtered to convex hull and time frame of species of interest

# Author: Naroa Olivia Schweizer
# last update: 04.01.26

# set correct working directory: setwd("C:/Users/Admin/Documents/Uni_Zuerich/Masters_Thesis/CoralReef_fish_local")

# load libraries
library(dplyr)          # data wrangling (filter, summarize, join)
library(stringr)        # consistent, readable string manipulation
library(sf)             # used to calculate geographic locations/convex hull
library(ggplot2)        # general-purpose plotting
library(missForest)     # randomforest to fill in NA values for env. data

#-----------------------------------------------------------------------------
## define paths and load files
#-----------------------------------------------------------------------------

## RECIFs
# define path to meta data for recifs
meta_recifs_path <- "intermed_outputs/RECIFs_inputs/meta_variables_recifs.csv"

# read meta data recifs
meta_recifs <- read.csv(meta_recifs_path)

# # define path to RECIFs zip file -> all data inputs
# data_recifs_path <- "intermed_outputs/RECIFs_inputs/RecifsDB.zip"
# 
# # only needed if not done in the first run already
# # unzip and load RECIFs data into a temp dir stored in R (will be deleted once session is closed)
# recifs_tmp <- tempdir()
# unzip(zipfile = data_recifs_path, exdir = recifs_tmp)

# otherwise just load: RECIFs_10kmbuffer_yearly_averages
# define path to REICFs yearly averages 10 km buffer
recifs_10km_yavg_path <- "intermed_outputs/RECIFs_inputs/recifs_data_010_yearlyaverag.RData"

# load R object
load(recifs_10km_yavg_path)

## Sample Info (to filter for geographic location)
# define path to run table (corrected for coherent coordinates!)
RunTable_path <-"data/EpiStri_RunTable.csv"

# read run table
RunTable <- read.csv(RunTable_path, header = TRUE, stringsAsFactors = FALSE)

## RECIFs structure: 
# 61'038 (5 x 5 km) reef raster cells (representing all reefs worldwide) 
# on reef -> buffer 2.5 km -> then also buffer for 10, 25 and 50 km (number at the end _002, _010 etc.)
# for each year: monthly average values

# create directory for adjusted env. data
dir.create('intermed_outputs/epinephelus_striatus/env_data_inputs', showWarnings = FALSE)

# define output path to save preped env data
save_path <- "intermed_outputs/epinephelus_striatus/env_data_inputs"

## define species specific end year for animation calculations
# define Yend
yend <- 2014

#=============================================================================
## 1. first output -> 10 km buffer, averages for each year
#=============================================================================

## only needs to be done with the OG files in the temp dir

## structure the RECIFs data
# make a list of all files from the temporary unzipped data
recifs_files <- list.files(recifs_tmp, recursive = TRUE, full.names = TRUE)

# exclude the png from the list and only include the .rda files
recifs_rda_files <- recifs_files[grepl("\\.rda$", recifs_files)]

# load all .rda files from RECIFs into one data list -> base data file
recifs_data <- setNames(
  lapply(recifs_rda_files, function(f) {
    e <- new.env()
    load(f, envir = e)
    e[[ls(e)[1]]]  # extract first object from each file (there is only one object in each .rda file)
  }),
  tools::file_path_sans_ext(basename(recifs_rda_files))
)

#-----------------------------------------------------------------------------
## 1.1) calculate yearly average out of monthly values
#-----------------------------------------------------------------------------

## spatial resolution of environmental data -> choose 10 km
recifs_data_010 <- recifs_data[
  grepl("_010$", names(recifs_data)) | grepl("coords", names(recifs_data))
]

## adjust the monthly values to a yearly average
# list the variable that have monthly values
vars_to_average <- c("CHL", "DHW", "FE", "O2", "PH", "NO3", "PO4", "SCV", "SPM", "SSS", "SST")

# prepare output list for yearly data
recifs_data_010_yearly <- list()

# loop through each variable product
for (v in vars_to_average) {
  
  obj_name <- paste0("s_allDB_", v, "_010")
  if (!obj_name %in% names(recifs_data_010)) next
  
  df <- as.data.frame(recifs_data_010[[obj_name]])
  
  # match columns ending with "_XMM.YYYY" or "_XMM_YYYY"
  data_cols <- grep(paste0(v, ".*_X[0-9]{2}[._][0-9]{4}$"), colnames(df), value = TRUE)
  if (length(data_cols) == 0) next
  
  # extract years
  years <- str_extract(data_cols, "[0-9]{4}$")
  
  # get unique years
  unique_years <- unique(years)
  
  # prepare new dataframe for yearly means
  yearly_df <- df[, grep("lon|lat", names(df), value = TRUE), drop = FALSE]
  
  # compute row-wise mean for each year (to not loose coordinates for each row when combining it later -> keeps row order the same)
  for (yr in unique_years) {
    cols_this_year <- data_cols[years == yr]
    yearly_df[[paste0(v, "_", yr)]] <- rowMeans(df[, cols_this_year, drop = FALSE], na.rm = TRUE)
  }
  
  recifs_data_010_yearly[[obj_name]] <- yearly_df
}

#-----------------------------------------------------------------------------
## 1.2) combine the averaged data list with the other variables in the list
#-----------------------------------------------------------------------------

# define static items
static_items <- c("allDB_coords",
                  "s_allDB_BATHY_010",
                  "s_allDB_LAND_010")
# define prefix for non static items, that don't need to be averaged
prefix_items <- c("s_allDB_CROP", "s_allDB_PDEN", "s_allDB_URBA", "s_allDB_VBD")

# add static items
for (item in static_items) {
  if (item %in% names(recifs_data_010)) {
    recifs_data_010_yearly[[item]] <- as.data.frame(recifs_data_010[[item]])
  }
}

# add all columns starting with a prefix
for (prefix in prefix_items) {
  cols <- grep(paste0("^", prefix), names(recifs_data_010), value = TRUE)
  for (col in cols) {
    recifs_data_010_yearly[[col]] <- as.data.frame(recifs_data_010[[col]])
  }
}

#-----------------------------------------------------------------------------
## 1.3) save the first output data as an r object 
#-----------------------------------------------------------------------------

save(recifs_data_010_yearly,
     file = file.path(save_path, "recifs_data_010_yearlyaverag.RData")
)

#=============================================================================
## 2) second output -> avg and sd for each time series variable
#=============================================================================

# this need do be implemented for each species separately

## filter the r object -> calculate mean and sd per time series variable
# list the variables with big time frames -> trim them to 2011 (for pterois volitans)
variable_info <- list(
  CHL = list(years = 1997:yend),
  DHW = list(years = 1985:yend),
  FE  = list(years = 1993:yend),
  O2  = list(years = 1993:yend),
  PH  = list(years = 1993:yend),
  NO3 = list(years = 1993:yend),
  PO4 = list(years = 1993:yend),
  SCV = list(years = 1993:yend),
  SPM = list(years = 1997:yend),
  SSS = list(years = 1993:yend),
  SST = list(years = 1985:yend)
)

# list other variables with short time frames
other_varCodes <- c("CROP", "PDEN", "URBA", "VBD")

# create result list
RECIFs_list <- list()

# loop through big time frame data and calculate avg and sd
for (v in names(variable_info)) {
  
  obj_name <- paste0("s_allDB_", v, "_010")
  
  # skip if variable does not exist
  if (!obj_name %in% names(recifs_data_010_yearly)) {
    warning("Missing yearly data object for ", v)
    next
  }
  
  df <- recifs_data_010_yearly[[obj_name]]
  
  # find columns like: CHL_1997, CHL_1998, …
  pattern <- paste0("^", v, "_[0-9]{4}$")
  matching_cols <- grep(pattern, colnames(df), value = TRUE)
  
  if (length(matching_cols) == 0) {
    warning("No yearly columns found for ", v)
    next
  }
  
  # extract the years
  years_found <- as.numeric(sub(".*_", "", matching_cols))
  
  # restrict to the time frame provided in variable_info
  wanted_years <- variable_info[[v]]$years
  use_cols <- matching_cols[years_found %in% wanted_years]
  
  if (length(use_cols) == 0) {
    warning("Columns exist for ", v, " but none match wanted years.")
    next
  }
  
  # compute stats
  df_subset <- df[, use_cols, drop = FALSE]
  
  avg <- rowMeans(df_subset, na.rm = TRUE)
  sd  <- apply(df_subset, 1, sd, na.rm = TRUE)
  
  # store into final list
  RECIFs_list[[paste0("s_allDB_", v, "_avg")]] <- avg
  RECIFs_list[[paste0("s_allDB_", v, "_sd")]]  <- sd
}

# loop through short time series to calculate avg & sd
for (v in other_varCodes) {
  
  matching_names <- grep(
    paste0("^s_allDB_", v),
    names(recifs_data_010_yearly),
    value = TRUE
  )
  
  if (length(matching_names) == 0) next
  
  df_combined <- do.call(
    cbind,
    lapply(recifs_data_010_yearly[matching_names], as.data.frame)
  )
  
  avg <- rowMeans(df_combined, na.rm = TRUE)
  sd  <- apply(df_combined, 1, sd, na.rm = TRUE)
  
  RECIFs_list[[paste0("s_allDB_", v, "_avg")]] <- avg
  RECIFs_list[[paste0("s_allDB_", v, "_sd")]]  <- sd
}

# define static items
static_items <- c("allDB_coords",
                  "s_allDB_BATHY_010",
                  "s_allDB_LAND_010")

# again add static items (same as above 1.2)
for (item in static_items) {
  if (item %in% names(recifs_data_010_yearly)) {
    RECIFs_list[[item]] <- as.data.frame(recifs_data_010_yearly[[item]])
  }
}

#-----------------------------------------------------------------------------
## 2.1) save the second output data as an r object 
#-----------------------------------------------------------------------------

# save the filtered data as an r object
save(RECIFs_list,
     file = file.path(save_path, "recifs_epistri_avg_sd.RData")
)

#=============================================================================
## 3) third output -> species specific csv of second output
#=============================================================================

#-----------------------------------------------------------------------------
## 3.1) define geographic boundaries
#-----------------------------------------------------------------------------

# Note: make sure the input Coordinate Reference System (CRS) is correct!

# convert sampling coordinates to geometry sf points add sample id as rownames
samp_sf <- st_as_sf(
  RunTable[, c("longitude", "latitude")],
  coords = c("longitude", "latitude"),
  crs = 4326
)
rownames(samp_sf) <- RunTable$Run

# compute convex hull polygon
convex_hull <- samp_sf %>%
  st_union() %>%      
  st_convex_hull()     

# convert to df to check if correctly computed
hull_df <- st_as_sf(convex_hull)
# plot convexhull
ggplot() +
  geom_sf(data = convex_hull, fill = NA, color = "grey60", size = 1) +
  geom_sf(data = samp_sf, color = "black", size = 2) +
  labs(
    x = "Longitude",
    y = "Latitude",
    title = "Convex Hull of Sampling Sites for Pterois volitans"
  ) +
  theme_minimal()

# check if pattern is the same as in the original paper!

#-----------------------------------------------------------------------------
## 3.2) define time frame
#-----------------------------------------------------------------------------

# define time frame
start_year <- 1980
end_year   <- yend

#-----------------------------------------------------------------------------
## 3.3) filter RECIFs data accordingly (time frame and geographic boundary)
#-----------------------------------------------------------------------------

# Note: here only the first output (O1) needs to be filtered for geographic 
# boundaries and time frame. This output will be used for animations to 
# facilitate interpretation.

## load r-object if already calculated
# # define path recifs_data_010_yearly (r-object)
# RECIFs_anim_path <- "intermed_outputs/RECIFs_inputs/recifs_data_010_yearlyaverag.RData"
# 
# # load r-object
# load(RECIFs_anim_path)

# O1: change RECIFs_list (the one with monthly averages) list into df
RECIFs_anim <- do.call(data.frame, recifs_data_010_yearly)

# convert RECIFs_anim to sf points (geometry)
RECIFs_anim_sf <- st_as_sf(
  RECIFs_anim,
  coords = c("allDB_coords.lon", "allDB_coords.lat"),
  crs = 4326,
  remove = FALSE
)

# O1: filter recifs rows to keep coordinates inside convex hull
anim_inside_mask <- st_within(RECIFs_anim_sf, convex_hull, sparse = FALSE)
anim_recifs_inside <- RECIFs_anim_sf[anim_inside_mask, ]

# O1: convert RECIFs_inside back to matrix format
RECIFs_anim <- anim_recifs_inside %>%
  st_drop_geometry()

# extract years from column names (if present)
col_names <- colnames(RECIFs_anim)
col_years <- suppressWarnings(
  as.integer(sub(".*?(\\d{4})$", "\\1", col_names))
)

# keep: columns with no year suffix and for those with -> the ones within time frame
keep_cols <- is.na(col_years) |
  (col_years >= start_year & col_years <= end_year)

RECIFs_anim <- RECIFs_anim[, keep_cols]

#-----------------------------------------------------------------------------
## 3.3) extract env. data points closest to sample sites (Nearest Neighbor)
#-----------------------------------------------------------------------------

# O2: change RECIFs_list (sd and averages for specific time frame) list into df
RECIFs_matrix <- do.call(data.frame, RECIFs_list)

# O2: convert RECIFs_matrix to sf points (geometry)
RECIFs_sf <- st_as_sf(
  RECIFs_matrix,
  coords = c("allDB_coords.lon", "allDB_coords.lat"),
  crs = 4326,
  remove = FALSE
)

# for each sample, find the nearest reef
nearest_idx <- st_nearest_feature(samp_sf, RECIFs_sf)

# extract the nearest reefs
nearest_reefs <- RECIFs_sf[nearest_idx, ] %>% st_drop_geometry()

# join the sample info (including sample name, lat, lon)
RECIFs_sample <- RunTable %>%
  select(Run, latitude, longitude) %>%
  bind_cols(nearest_reefs) %>%
  rename(
    sample_lon = longitude,
    sample_lat = latitude,
    sample = Run
  )

#-----------------------------------------------------------------------------
## 3.4) replace NA values in environmental matrix using missForest function
#-----------------------------------------------------------------------------

# make column sample the row names (same as for imputed_gt)
rownames(RECIFs_sample) <- RECIFs_sample$sample

# exclude sample as column (as this has charater and input cells need to be num)
RECIFs_sample <- RECIFs_sample %>%
  select(-sample)
# check if all are numeric using -> sapply(RECIFs_sample, class)

# replace NA values using random forest algorithm 
set.seed(42) # for reproducibility set seed point

imputed_env <- missForest(
  RECIFs_sample,
  maxiter = 10,
  ntree = 300,
  parallelize = "no"
)

RECIFs_sample <- imputed_env$ximp

#-----------------------------------------------------------------------------
## 3.5) extract env. data points closest to sample sites (Nearest Neighbor)
#-----------------------------------------------------------------------------

# ## load r-object if already calculated (or loaded)
# # define path recifs_ptevol_avg_sd.RData (r-object)
# RECIFs_list_path <- "intermed_outputs/siphamia_tubifer/env_data_inputs/recifs_siphtub_avg_sd.RData"
# 
# # load r-object
# load(RECIFs_list_path)

# change RECIFs_list (sd and averages for specific time frame) list into df
RECIFs_matrix <- do.call(data.frame, RECIFs_list)

# convert RECIFs_matrix to sf points (geometry)
RECIFs_matrix_sf <- st_as_sf(
  RECIFs_matrix,
  coords = c("allDB_coords.lon", "allDB_coords.lat"),
  crs = 4326,
  remove = FALSE
)

# filter recifs rows to keep coordinates inside convex hull
matrix_inside_mask <- st_within(RECIFs_matrix_sf, convex_hull, sparse = FALSE)
matrix_recifs_inside <- RECIFs_matrix_sf[matrix_inside_mask, ]

# convert RECIFs_inside back to matrix format
RECIFs_matrix <- matrix_recifs_inside %>%
  st_drop_geometry()

# replace NA values using random forest algorithm 
set.seed(100) # for reproducibility set seed point

imputed_matrix <- missForest(
  RECIFs_matrix,
  maxiter = 10,
  ntree = 300,
  parallelize = "no"
)

# replace imputed values back into RECIFs_matrix
RECIFs_matrix <- imputed_matrix$ximp

#-----------------------------------------------------------------------------
## 3.6) save different environmental data outputs
#-----------------------------------------------------------------------------

# save filtered RECIFs_anim data -> monthly averages per year within convex hull
# and time frame
write.csv(RECIFs_anim, file = file.path(save_path, "EpiStri_RECIFs_anim.csv"), 
          row.names = FALSE)

# save filtered RECIFs_sample data -> overall sd and avg of time frame filtered
# for reef cell closest (NN) to sample coordinates, imputed for NA values using 
# missForest
write.csv(RECIFs_sample, file = file.path(save_path, "EpiStri_RECIFs_sample.csv"), 
          row.names = FALSE)

# save filtered RECIFs_matrix data -> imputed NA values of reef cells using 
# missForest function
write.csv(RECIFs_matrix, file = file.path(save_path, "EpiStri_RECIFs_matrix.csv"), 
          row.names = FALSE)
