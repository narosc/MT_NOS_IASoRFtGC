#=============================================================================
### first analysis - prepare genomic data for analysis
#=============================================================================

# This script prepares the genomic data. It imputes the missing data using LEA 
# package snmf function. Futhermore, it excludes the individuals with no 
# geolocation and those with insufficient quality reads.
# Inputs:
# 1) compressed vcf file (already filtered for MAF, missing GTs and biallelic)
# 2) run table with metadata for samples downloaded from NCBI Bioproject
# Outputs: 
# 1. a filtered genlight object only including high quality samples with geolocation
# 2. an adjusted/filtered run table for species

# Author: Naroa Olivia Schweizer
# last update: 02.03.26

# set correct working directory: setwd("C:/Users/Admin/Documents/Uni_Zuerich/Masters_Thesis/CoralReef_fish_local")

# load libraries
library(dplyr)          # data wrangling (filter, summarize, join)
library(vcfR)           # read and manipulate VCF files
library(adegenet)       # PCA & multivariate analysis of genetic data

#-----------------------------------------------------------------------------
## define paths and load files
#-----------------------------------------------------------------------------

## Sample Info (to check if geolocation is included)
# define path to run table
RunTable_path <-"data/SraRunTable_PRJEB36904_Ephinephelusstriatus.csv"

# read run table
RunTable <- read.csv(RunTable_path, header = TRUE, stringsAsFactors = FALSE)

## filtered VCF file (biallelic, MAF > 0.05, sites with less than 10% missing genotypes)
# define path to filtered vcf
vcf_path <- "intermed_outputs/epinephelus_striatus/var/cohort_biallelic_fmiss_maf_pruned.vcf.gz"

# read VCF
vcf <- vcfR::read.vcfR(vcf_path)

## outputs
# define output path
output_path1 <- "intermed_outputs/epinephelus_striatus/var"
output_path2 <- "data"

#=============================================================================
## 1) prepare and filter variants and RunTable
#=============================================================================

#-----------------------------------------------------------------------------
## 1.1) convert compressed filtered vcf into genlight object  
#-----------------------------------------------------------------------------

# convert VCF to a genlight object
genlight_obj <- vcfR::vcfR2genlight(vcf)

# clean genlight sample names
inds <- indNames(genlight_obj)              # get sample names from genlight
inds <- basename(inds)                      # strip path
inds <- sub("\\.sorted\\.bam$", "", inds)   # strip extension
indNames(genlight_obj) <- inds              # assign back

#-----------------------------------------------------------------------------
## 1.2) filter RunTable for quality-sufficient individuals with geolocations
#-----------------------------------------------------------------------------

# Note: for Pterois volitans individual SRR8354500 and SRR8354521 have 
# insufficient quality (weird GC content). Filter them out.
# no need for second species Amphiprion bicinctus!!
# for Suphaima tubifer read SRR5489654 has a weird GC content, filter out.

# exclude bad quality individuals
RunTable <- RunTable %>%
  filter(!Run %in% c("ERR3950796", "ERR3950797", "ERR3950798",
                     "ERR3950829", "ERR3950836", "ERR3950864"))

# notes for Siphamia tubifer:
# RunTable <- RunTable %>%
#   filter(!Run %in% c("SRR5489527", "SRR5489521", "SRR5489526", "SRR5489528", "SRR5489654"))
# # kinda went missing in alignment -> 2 in Yonabaru, 2 in Sesoko -> abnormal GC content ("SRR5489654")
# # notes for Pterois volitans
# RunTable <- RunTable %>%
#   filter(!Run %in% c("SRR8354500", "SRR8354521"))
# # notes for Epinephelus striatus
# RunTable <- RunTable %>%
#   filter(!Run %in% c("ERR3950796", "ERR3950797", "ERR3950798",
#                     "ERR3950829", "ERR3950836", "ERR3950864"))
# either went missing in computation or no mapped reads (0% rate) -> "ERR3950829".

# exclude individuals with missing geolocation
RunTable <- RunTable %>%
  filter(!is.na(latitude) & !is.na(longitude))

#-----------------------------------------------------------------------------
## 1.3) filter genlight object individuals (sufficient quality, geolocation)
#-----------------------------------------------------------------------------

# individuals to keep (from RunTable)
keep_inds <- RunTable$Run

# filter genlight object to only keep_inds
genlight_obj <- genlight_obj[indNames(genlight_obj) %in% keep_inds, ]

#=============================================================================
## 2) save outputs
#=============================================================================

# save filtered genlight object
saveRDS(genlight_obj,
        file = file.path(output_path1, "estriatus_genlight_filtered.rds"))

# save filtered RunTable
write.csv(RunTable, file = file.path(output_path2, "EpiStri_RunTable.csv"), 
          row.names = FALSE)
