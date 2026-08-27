#=============================================================================
### eighth analysis - Protein to SNP - Preparation for PicMin
#=============================================================================

# This script maps the SNPs to the protein annotation of a reference genome 
# (the output from the BLAST -> zebra fish proteins associated to reference genome).
# Inputs: 
# 1) SNPs list (R Object) per environmental variable (output: LFMM GEA association)
# 2) BLAST protein associated reference genome
# Output: 
# 1. R Object: List with protein_id, combined p-value (using stouffer method), 
# number of snps on protein and their "name" per environmental variable (serving
# as the input for PicMin)

# Author: Naroa Olivia Schweizer
# last update: 11.06.26

# set correct working directory: setwd("C:/Users/Admin/Documents/Uni_Zuerich/Masters_Thesis/CoralReef_fish_local")

# load libraries
library(dplyr)          # data wrangling (filter, summarize, join)
library(stringr)        # consistent, readable string manipulation
library(poolr)          # use stouffer method to get one p-value per protein
library(tidyverse)      # restructure tables and formats

#-----------------------------------------------------------------------------
## define paths and load files
#-----------------------------------------------------------------------------

## R-Object of SNPs per env. variable LFMM
# define path to r-object
SNPs_list_path <- "intermed_outputs/epinephelus_striatus/lfmm/lfmm_results.RData"

# load r-object (results_list)
load(SNPs_list_path)

## R-Object of SNPs per env. variable CA
# define path to r-object
SNPs_list_path <- "intermed_outputs/epinephelus_striatus/ca/ca_results.RData"

# load r-object (results_list)
load(SNPs_list_path)

## protein associated reference genome
# define path to protein ref genome
prot_ref_path <- "intermed_outputs/epinephelus_striatus/protein_annotation/drerio_vs_estriatus.tblastn.out"

# load prot ref
prot_ref <- read.table(prot_ref_path, sep = "\t", header = FALSE)

# add column names for the prot_ref (defined in slurm script)
colnames(prot_ref) <- c(
  "protein_id",
  "contig",
  "percent_identity",
  "alignment_length",
  "mismatches",
  "gap_opens",
  "query_start",
  "query_end",
  "genome_start",
  "genome_end",
  "evalue",
  "bitscore"
)

## Danio rerio protein id and gene id information text file
# define path to gene id and protein id file
gtf_zip <- "intermed_outputs/danio_rerio/Danio_rerio.GRCz11.115.chr.gtf.gz"

# create a temporary directory
tmp_dir <- tempdir()

# unzip into the temporary directory
gtf_path <- R.utils::gunzip(
  filename = gtf_zip,
  destname = file.path(tmp_dir, "Danio_rerio.GRCz11.115.chr.gtf"),
  overwrite = TRUE,
  remove = FALSE
)

# load prot ref
gtf <- read.delim(
  gtf_path,
  header = FALSE,
  sep = "\t",
  comment.char = "#",
  quote = "",
  stringsAsFactors = FALSE
)

# add column names to the gtf
colnames(gtf) <- c(
  "seqname", "source", "feature", "start", "end",
  "score", "strand", "frame", "attribute"
)

# create df with only protein id, gene id, transcript id & gene name (note: information of protein to gene only within "CDS")
prot_gene_df <- gtf %>%
  filter(feature == "CDS") %>%
  mutate(
    gene_id = str_match(attribute, 'gene_id "([^"]+)"')[,2],
    transcript_id = str_match(attribute, 'transcript_id "([^"]+)"')[,2],
    protein_id = str_match(attribute, 'protein_id "([^"]+)"')[,2],
    gene_name = str_match(attribute, 'gene_name "([^"]+)"')[,2]
  ) %>%
  filter(!is.na(protein_id)) %>%
  select(gene_id, gene_name, transcript_id, protein_id) %>%
  distinct()

#=============================================================================
## 1) identify which SNPs sit on which protein, reconstruct df
#=============================================================================

# filter out artefact snps for Epinephelus striatus -> 22 SNPs have NA (NaN) 
# values in both LFMM and CA for all environmental variables -> all are 
# heterozygous SNPs, does don't show any variation

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

# remove monomorphic heterozygous SNPs from all environmental variables
results_list <- lapply(results_list, function(df) {
  
  df <- df[!df$SNP %in% NA_NaN_snps, ]
  
  return(df)
})

#-----------------------------------------------------------------------------
## 1.1) prepare protein ref 
#-----------------------------------------------------------------------------

# make sure start and end in bp are correct -> add new columns with start end bp
prot_ref <- prot_ref %>%
  mutate(
    g_start = pmin(genome_start, genome_end),
    g_end   = pmax(genome_start, genome_end)
  )

# filter protein reference to only include high quality hits (avoid spurious hits from blast output)
prot_ref_fil <- prot_ref %>%
  filter(
    percent_identity > 70
  )

#-----------------------------------------------------------------------------
## 1.2) extract SNPs position from contig - function
#-----------------------------------------------------------------------------

# # function to separate SNP-chromosome and position (this worked for Amphiprion bicinctus -> because not contig but chromosome)
# extract_snp_pos <- function(df){
#   df %>%
#     mutate(
#       contig = str_extract(SNP, "^.+(?=_[^_]+$)"),
#       pos    = as.numeric(str_extract(SNP, "[^_]+$"))
#     )
# }

# function to separate SNP-contig and position (this worked for Siphamia tubifer, Pterois volitans, Epinephelus striatus)
extract_snp_pos <- function(df){
  df %>%
    mutate(
      contig = str_extract(SNP, "^[^_]+"),
      pos    = as.numeric(str_extract(SNP, "[0-9]+$"))
    )
}

#-----------------------------------------------------------------------------
## 1.3) map SNPs to protein and then retreive gene_id(s)
#-----------------------------------------------------------------------------

# function matching snps to protein (window +-10'000 bp) and then retreive gene_id
join_snps_to_genes <- function(df, prot_ref_fil, prot_gene_df) {
  # create empty gene column first
  df$gene <- ""
  # set rownames outside of the loop
  rownames(prot_gene_df) <- prot_gene_df$protein_id
  # loop over rows
  for (i in seq_len(nrow(df))) {
    contig <- df$contig[i]
    pos <- df$pos[i]
    # subset proteins near SNP (linkage disequillibrium decay considered)
    s_prot_ref <- prot_ref_fil[
      prot_ref_fil$contig == contig &
        (prot_ref_fil$g_start - 10000) <= pos &
        (prot_ref_fil$g_end + 10000) >= pos,
    ]
    # skip if no matches (no nearby protein)
    if (nrow(s_prot_ref) == 0) {
      next
    }
    # remove version suffix (last two characters)
    s_prot_ref$protein_id2 <- substr(s_prot_ref$protein_id, 1, nchar(s_prot_ref$protein_id) - 2)
    # map protein IDs to gene IDs
    prots <- prot_gene_df[s_prot_ref$protein_id2, ]
    # create output
    output <- paste(
      unique(na.omit(prots$gene_id)),
      collapse = ","
    )
    # assign result to SAME ROW
    df$gene[i] <- output
  }
  return(df)
}

#-----------------------------------------------------------------------------
## 1.4) combine p-values of snps per gene
#-----------------------------------------------------------------------------

# Note: stouffer method converts p-values to normal z-scores, then averages the
# z-scores and then converts this averaged z-score back to a p-value...

# 1. duplicate snps rows so that  each gene in the gene column is seperate 
# 2. filter out emtpy strings
# 3. group by gene (gene_id)
# 4. calculate combined p-value for each group (do not collapse any columns) 
# -> so that each gene has a unique combined p-value

# LFMM output version
combine_gene_pvals_LFMM <- function(df) {
  # split rows by gene
  df_sep <- df %>%
    separate_rows(gene, sep = ",")
  # calculate combined p-value per gene (exclude rows with emtpy string -> no genes assigned)
  gene_pvals <- df_sep %>%
    filter(gene != "") %>%
    group_by(gene) %>%
    summarise(
      combined_p = stouffer(p_gif)$p,
      .groups = "drop"
    )
  # join combined p-values back to gene seperated df
  df_final <- df_sep %>%
    left_join(gene_pvals, by = "gene")

  return(df_final)
}

# CA output version
combine_gene_pvals_CA <- function(df) {
  # split rows by gene
  df_sep <- df %>%
    separate_rows(gene, sep = ",")
  # calculate combined p-value per gene (exclude rows with emtpy string -> no genes assigned)
  gene_pvals <- df_sep %>%
    filter(gene != "") %>%
    group_by(gene) %>%
    summarise(
      combined_p = stouffer(pvalue)$p,
      .groups = "drop"
    )
  # join combined p-values back to gene seperated df
  df_final <- df_sep %>%
    left_join(gene_pvals, by = "gene")

  return(df_final)
}

#-----------------------------------------------------------------------------
## 1.5) wrap steps and loop over all env. variables separately - function
#-----------------------------------------------------------------------------

## LFMM version
# wrap steps from above
process_env_variable <- function(df, prot_ref_fil, prot_gene_df){
  # extract contig and SNP position
  df_pos <- extract_snp_pos(df)
  # map SNPs to proteins and then to genes
  mapped <- join_snps_to_genes(
    df_pos,
    prot_ref_fil,
    prot_gene_df
  )
  # calculate combined gene-level p-values
  gene_table <- combine_gene_pvals_LFMM(mapped)
  return(gene_table)
}

## CA version
# wrap steps from above
process_env_variable <- function(df, prot_ref_fil, prot_gene_df){
  # extract contig and SNP position
  df_pos <- extract_snp_pos(df)
  # map SNPs to proteins and then to genes
  mapped <- join_snps_to_genes(
    df_pos,
    prot_ref_fil,
    prot_gene_df
  )
  # calculate combined gene-level p-values
  gene_table <- combine_gene_pvals_CA(mapped)
  return(gene_table)
}

# loop over results_list to get an r-object of 1 p-value per row
# and a gene-level gif-corrected p-value per column
gene_results <- lapply(
  results_list,
  process_env_variable,
  prot_ref_fil = prot_ref_fil,
  prot_gene_df = prot_gene_df
)

#=============================================================================
## 2) save output r-object -> protein_list
#=============================================================================

## protein mapping results
# # LFMM
# # save SNPs mapping to genes results as r-object
# save(gene_results, file = "intermed_outputs/epinephelus_striatus/protein_annotation/SNPs_mapped_to_genes_LFMM.RData")
# CA
# save SNPs mapping to genes results as r-object
save(gene_results, file = "intermed_outputs/epinephelus_striatus/protein_annotation/SNPs_mapped_to_genes_CA.RData")
