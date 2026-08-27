## notes: the input tables are not stored on the cluster, can just download them 
# from NCBI... have them stored on google drive (save storage)

# Author: Naroa Olivia Schweizer
# last update: 29.04.26

# load libraries
library(dplyr)
library(tidyr)
library(stringr)
library(tibble)

#-----------------------------------------------------------------------------
## Ephinephelus striatus
#-----------------------------------------------------------------------------

## what has been done for Ephinephelus striatus
# define path to run table
NoGeoLoc_path <-"data/SraRunTable_PRJEB36904_Ephinephelusstriatus_OG.csv"

# read run table
NoGeoLoc <- read.csv(NoGeoLoc_path, header = TRUE, stringsAsFactors = FALSE)

# geolocations were retrieved like so: google maps > name of the reef 
# (geographic_location_(region_and_locality)) > coordinates for approximate 
# position on map in paper is retrieved and defined here
geolocations_bahamas <- data.frame(
  reef = c(
    "Abaco", "Andros", "Eleuthera", "Exuma",
    "Great Inagua", "Hail Mary", "Long Island",
    "New Providence", "Ragged Island"
  ),
  latitude = c(
    26.296740, 24.528850, 25.213777, 24.013479,
    21.125450, 23.581889, 23.114118,
    25.032612, 22.399288
  ),
  longitude = c(
    -77.534347, -77.701774, -76.119112, -76.246584,
    -73.667139, -75.353563, -74.945766,
    -77.564033, -75.792046
  ),
  stringsAsFactors = FALSE
)

# join tables
GeoLocBah <- NoGeoLoc %>%
  left_join(
    geolocations_bahamas,
    by = c("geographic_location_.region_and_locality." = "reef")
  )

# save filtered RunTable
write.csv(GeoLocBah, file = file.path("data/SraRunTable_PRJEB36904_Ephinephelusstriatus.csv"), 
          row.names = FALSE)

#-----------------------------------------------------------------------------
## Dascylltus trimaculatus
#-----------------------------------------------------------------------------

## what has been done for Dascylltus trimaculatus
# define path to run table
NoGeoLoc_path <-"data/SraRunTable_PRJNA473382_Dascyllustrimaculatus_OG.csv"

# read run table
NoGeoLoc <- read.csv(NoGeoLoc_path, header = TRUE, stringsAsFactors = FALSE)

# separate lat and long into a column each
WithGeoLoc <- NoGeoLoc %>%
  separate(lat_lon, into = c("lat", "lat_dir", "lon", "lon_dir"), sep = " ") %>%
  mutate(
    latitude = as.numeric(lat) * if_else(lat_dir == "S", -1, 1),
    longitude = as.numeric(lon) * if_else(lon_dir == "W", -1, 1)
  ) %>%
  select(-lat_dir, -lon_dir, -lat, -lon)

# save filtered RunTable
write.csv(WithGeoLoc, file = file.path("data/SraRunTable_PRJNA473382_Dascyllustrimaculatus.csv"), 
          row.names = FALSE)

#-----------------------------------------------------------------------------
## Siphamia tubifer
#-----------------------------------------------------------------------------

## what has been done for Siphamia tubifer
# define path to run table
NoGeoLoc_path <-"data/SraRunTable_PRJNA385011_Siphamiatubifer_OG_noGEOloc.csv"

# read run table
NoGeoLoc <- read.csv(NoGeoLoc_path, header = TRUE, stringsAsFactors = FALSE)

# define path to join geo loc
GeoLoc_path <-"data/PRJNA385011_Siphamiatubifer_GeoLocations.csv"

# read geo loc table
GeoLoc <- read.csv(GeoLoc_path, header = TRUE, stringsAsFactors = FALSE)

# add new column for the join with the stud_names_prefix to run table
NoGeoLoc <- NoGeoLoc %>%
  mutate(study_book_prefix = sub("[0-9].*", "", stud_book_number))

# left join on tables to add geographic information
merged <- NoGeoLoc %>%
  left_join(GeoLoc, by = c("study_book_prefix" = "stud_book_prefix"))

# save filtered RunTable
write.csv(merged, file = file.path("data/SraRunTable_PRJNA385011_Siphamiatubifer.csv"), 
          row.names = FALSE)
