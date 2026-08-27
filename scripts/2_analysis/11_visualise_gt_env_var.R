#==============================================================================
### eleventh analysis - geographically visualise PicMin outpus
#==============================================================================

# This script visualises significant GTs by mapping them onto the heatmap of
# the associated env variable
# Inputs:
# 1) genlight object of each species
# 2) run table of each species
# 3) environmental input (RECIFs_species files)
# 4) mapped SNPs to gene files for each species
# Output:
# 1. TWO summary figures (2 species each, alphabetical order), each showing a
#    map with per-site allele frequency (of whichever method's SNP had the
#    lower p-value) + env. variable, and two GT vs env. variable panels (one
#    for LFMM, one for CA, each with its own rounded p-value in the subtitle)
# 2. CSV of every gene (both methods, all env. variables) with pooled_q < 0.25,
#    with the p-value of the "most contributing" SNP per species

# set correct working directory:
setwd("C:/Users/Admin/Documents/Uni_Zuerich/Masters_Thesis/CoralReef_fish_local")

# Author: Naroa Olivia Schweizer
# last update: 21.08.26

# load libraries
library(dplyr)          # data wrangling (filter, summarize, join)
library(tibble)         # rownames_to_column etc.
library(ggplot2)        # general-purpose plotting
library(adegenet)       # genlight object methods (as.matrix on genlight)
library(ggspatial)      # scale bar / north arrow map annotations
library(grid)           # needed for unit() (colourbar sizing)
library(ggnewscale)     # for multiple colour scales in one plot
library(rnaturalearth)  # for clean world map polygons
library(purrr)          # iterating over lists/vectors without writing for loops
library(patchwork)      # combining multiple plots

#------------------------------------------------------------------------------
## define paths and load files
#------------------------------------------------------------------------------

# define paths in a list to all needed files for the 4 species
species_list <- list(
  ptevol = list(
    species             = "Pterois volitans",
    genlight_path       = "intermed_outputs/pterois_volitans/var/ptevol_genlight_filtered_imputed.rds",
    runtable_path       = "data/PteVol_RunTable.csv",
    recifs_matrix_path  = "intermed_outputs/pterois_volitans/env_data_inputs/PteVol_RECIFs_matrix.csv",
    recifs_sample_path  = "intermed_outputs/pterois_volitans/env_data_inputs/PteVol_RECIFs_sample.csv",
    lfmm_path           = "intermed_outputs/pterois_volitans/lfmm/lfmm_results.RData",
    ca_path             = "intermed_outputs/pterois_volitans/ca/ca_results.RData",
    gene_path_lfmm      = "intermed_outputs/pterois_volitans/protein_annotation/SNPs_mapped_to_genes_LFMM.RData",
    gene_path_ca        = "intermed_outputs/pterois_volitans/protein_annotation/SNPs_mapped_to_genes_CA.RData",
    xlim = c(-90, -60), ylim = c(10, 30)
  ),
  amphibic = list(
    species             = "Amphiprion bicinctus",
    genlight_path       = "intermed_outputs/amphiprion_bicinctus/var/amphibic_genlight_filtered_imputed.rds",
    runtable_path       = "data/AmphiBic_RunTable.csv",
    recifs_matrix_path  = "intermed_outputs/amphiprion_bicinctus/env_data_inputs/AmphiBic_RECIFs_matrix.csv",
    recifs_sample_path  = "intermed_outputs/amphiprion_bicinctus/env_data_inputs/AmphiBic_RECIFs_sample.csv",
    lfmm_path           = "intermed_outputs/amphiprion_bicinctus/lfmm/lfmm_results.RData",
    ca_path             = "intermed_outputs/amphiprion_bicinctus/ca/ca_results.RData",
    gene_path_lfmm      = "intermed_outputs/amphiprion_bicinctus/protein_annotation/SNPs_mapped_to_genes_LFMM.RData",
    gene_path_ca        = "intermed_outputs/amphiprion_bicinctus/protein_annotation/SNPs_mapped_to_genes_CA.RData",
    xlim = c(30, 60), ylim = c(10, 35)
  ),
  siphtub = list(
    species             = "Siphamia tubifer",
    genlight_path       = "intermed_outputs/siphamia_tubifer/var/siphtub_genlight_filtered_imputed.rds",
    runtable_path       = "data/SiphTub_RunTable.csv",
    recifs_matrix_path  = "intermed_outputs/siphamia_tubifer/env_data_inputs/SiphTub_RECIFs_matrix.csv",
    recifs_sample_path  = "intermed_outputs/siphamia_tubifer/env_data_inputs/SiphTub_RECIFs_sample.csv",
    lfmm_path           = "intermed_outputs/siphamia_tubifer/lfmm/lfmm_results.RData",
    ca_path             = "intermed_outputs/siphamia_tubifer/ca/ca_results.RData",
    gene_path_lfmm      = "intermed_outputs/siphamia_tubifer/protein_annotation/SNPs_mapped_to_genes_LFMM.RData",
    gene_path_ca        = "intermed_outputs/siphamia_tubifer/protein_annotation/SNPs_mapped_to_genes_CA.RData",
    xlim = c(126.75, 128.5), ylim = c(26, 27)
  ),
  epistri = list(
    species             = "Epinephelus striatus",
    genlight_path       = "intermed_outputs/epinephelus_striatus/var/epistri_genlight_filtered_imputed.rds",
    runtable_path       = "data/EpiStri_RunTable.csv",
    recifs_matrix_path  = "intermed_outputs/epinephelus_striatus/env_data_inputs/EpiStri_RECIFs_matrix.csv",
    recifs_sample_path  = "intermed_outputs/epinephelus_striatus/env_data_inputs/EpiStri_RECIFs_sample.csv",
    lfmm_path           = "intermed_outputs/epinephelus_striatus/lfmm/lfmm_results.RData",
    ca_path             = "intermed_outputs/epinephelus_striatus/ca/ca_results.RData",
    gene_path_lfmm      = "intermed_outputs/epinephelus_striatus/protein_annotation/SNPs_mapped_to_genes_LFMM.RData",
    gene_path_ca        = "intermed_outputs/epinephelus_striatus/protein_annotation/SNPs_mapped_to_genes_CA.RData",
    xlim = c(-79, -72), ylim = c(20.5, 27.5)
  )
)

# create function to load the species data list
load_species_data <- function(meta) {
  run_table <- read.csv(meta$runtable_path, header = TRUE, stringsAsFactors = FALSE)
  recifs_sample <- read.csv(meta$recifs_sample_path)
  rownames(recifs_sample) <- run_table$Run
  # important: both ca and lfmm results are called the same -> thus load in dumps
  # to not overwrite!
  lfmm_env <- new.env()
  load(meta$lfmm_path, envir = lfmm_env)
  ca_env <- new.env()
  load(meta$ca_path, envir = ca_env)
  gene_env_lfmm <- new.env()
  load(meta$gene_path_lfmm, envir = gene_env_lfmm)
  gene_env_ca <- new.env()
  load(meta$gene_path_ca, envir = gene_env_ca)
  # create structure of the list
  list(
    species       = meta$species,
    genlight      = readRDS(meta$genlight_path),
    run_table     = run_table,
    recifs_matrix = read.csv(meta$recifs_matrix_path),
    recifs_sample = recifs_sample,
    lfmm_results  = lfmm_env$results_list,
    ca_results    = ca_env$results_list,
    gene_lfmm     = gene_env_lfmm$gene_results,
    gene_ca       = gene_env_ca$gene_results,
    xlim          = meta$xlim,
    ylim          = meta$ylim
  )
}

# load all four species into one list, easy access for downstream analysis
all_data <- purrr::map(species_list, load_species_data)

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

# function to drop the 22 SNPs in the LFMM and CA results only for E. striatus
drop_snps <- function(results_list, snp_ids, id_col = "snp") {
  strip_one <- function(df) {
    if (!is.data.frame(df)) return(df)   # not a data.frame -> leave untouched, check manually
    if (id_col %in% names(df)) {
      df[!df[[id_col]] %in% snp_ids, , drop = FALSE]
    } else {
      df[!rownames(df) %in% snp_ids, , drop = FALSE]  # fall back to rownames as SNP id
    }
  }
  if (is.data.frame(results_list)) strip_one(results_list) else purrr::map(results_list, strip_one)
}

# apply function for E. striatus
all_data$epistri$lfmm_results <- drop_snps(all_data$epistri$lfmm_results, NA_NaN_snps)
all_data$epistri$ca_results   <- drop_snps(all_data$epistri$ca_results,   NA_NaN_snps)

## picmin results
# LFMM
picmin_path <- "intermed_outputs/picmin/picmin_results_lfmm.rda"
# load r-object (picmin_results)
load(picmin_path)
lfmm_picmin <- picmin_results

# CA
picmin_path <- "intermed_outputs/picmin/picmin_results_ca.rda"
# load r-object (picmin_results)
load(picmin_path)
ca_picmin <- picmin_results

# load world base map
world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")

#==============================================================================
## 1) find the "best gene" per environmental variable + its target SNP
#==============================================================================

# create list to store best picmin outputs (both LFMM and CA, and all species have same env. variables)
env_vars <- names(lfmm_picmin)

# function to retrieve the best gene -> q-value < 0.25
select_best_gene <- function(env_var, picmin_obj, method, q_threshold = 0.25) {
  df <- picmin_obj[[env_var]]
  df |>
    dplyr::filter(pooled_q < q_threshold) |>
    dplyr::arrange(pooled_q) |>                  # rank by strongest convergent signal
    dplyr::slice(1) |>                           # keep only the single best gene
    dplyr::transmute(locus, pooled_q, n_est, method = method)
}

# select "best gene" for LFMM -> lowest pooled_q < 0.25, any n_est
best_genes_lfmm <- purrr::map_dfr(env_vars, function(ev) {
  select_best_gene(ev, lfmm_picmin, method = "LFMM") |>
    dplyr::mutate(env_var = ev, .before = 1)
})

# select "best gene" for CA -> lowest pooled_q < 0.25, any n_est
best_genes_ca <- purrr::map_dfr(env_vars, function(ev) {
  select_best_gene(ev, ca_picmin, method = "CA") |>
    dplyr::mutate(env_var = ev, .before = 1)
})

# selected "best gene" will be searched within SNPs to genes file, to retrieve
# SNPs that sit on that gene, then in the GEA output looks for the SNP with the
# lowest q-value and selects that SNP as "most contributing" SNP
# find "most contributing" SNP on gene
best_snp_for_gene <- function(gene_id, env_var, gene_mapping_list, raw_results_list,
                              method = c("LFMM", "CA")) {
  method <- match.arg(method)
  rank_col <- switch(method, LFMM = "p_gif", CA = "pvalue")
  na_result <- list(snp = NA_character_, p = NA_real_)
  gene_df <- gene_mapping_list[[env_var]]
  if (is.null(gene_df)) return(na_result)
  candidate_snps <- gene_df$SNP[gene_df$gene == gene_id]
  if (length(candidate_snps) == 0) return(na_result)
  raw_df <- raw_results_list[[env_var]]
  if (is.null(raw_df)) return(na_result)
  hits <- raw_df[raw_df$SNP %in% candidate_snps, , drop = FALSE]
  if (nrow(hits) == 0) return(na_result)
  best_idx <- which.min(hits[[rank_col]])
  list(snp = hits$SNP[best_idx], p = hits[[rank_col]][best_idx])
}

# function to implement what's written above
get_target_snps <- function(gene_id, env_var, method, all_data,
                            gene_results_field, raw_results_field) {
  purrr::map_dfr(all_data, function(d) {
    res <- best_snp_for_gene(
      gene_id = gene_id, env_var = env_var,
      gene_mapping_list = d[[gene_results_field]],
      raw_results_list  = d[[raw_results_field]],
      method = method
    )
    tibble::tibble(
      species = d$species, env_var = env_var, method = method,
      gene = gene_id, snp = res$snp, p_value = res$p
    )
  })
}

#------------------------------------------------------------------------------
# 1.1) export CSV: every qualifying gene (q < 0.25), both methods, all env vars
#------------------------------------------------------------------------------

# function to collect all qualifying hits for one env var + method
collect_all_hits <- function(env_var, picmin_obj, method, q_threshold = 0.25) {
  df <- picmin_obj[[env_var]]
  df |>
    dplyr::filter(pooled_q < q_threshold) |>
    dplyr::transmute(locus, pooled_q, n_est, method = method)
}

# run for both methods, across all env vars
all_hits_lfmm <- purrr::map_dfr(env_vars, function(ev) {
  collect_all_hits(ev, lfmm_picmin, method = "LFMM") |>
    dplyr::mutate(env_var = ev, .before = 1)
})
all_hits_ca <- purrr::map_dfr(env_vars, function(ev) {
  collect_all_hits(ev, ca_picmin, method = "CA") |>
    dplyr::mutate(env_var = ev, .before = 1)
})
all_hits <- dplyr::bind_rows(all_hits_lfmm, all_hits_ca)

# for every qualifying (env_var, gene, method) hit, find the "most
# contributing" SNP
all_hits_snps <- purrr::pmap_dfr(all_hits, function(env_var, locus, pooled_q, n_est, method) {
  get_target_snps(
    gene_id = locus, env_var = env_var, method = method, all_data = all_data,
    gene_results_field = if (method == "LFMM") "gene_lfmm" else "gene_ca",
    raw_results_field  = if (method == "LFMM") "lfmm_results" else "ca_results"
  ) |>
    dplyr::mutate(pooled_q = pooled_q, n_est = n_est)
})

# final column order: env variable, gene id, pooled q-value, n_est, method, SNP, p-value
all_hits_snps <- all_hits_snps |>
  dplyr::select(env_var, gene, pooled_q, n_est, method, species, snp, p_value) |>
  dplyr::arrange(env_var, method, pooled_q, species)

# inspect
all_hits_snps

# save gene hit and p-value per species of most contributing SNP
write.csv(all_hits_snps, file = "outputs/picmin_4s/picmin_target_snps_q025.csv", row.names = FALSE)

#==============================================================================
## 2) per-site allele frequency / environmental value helpers
#==============================================================================

# Note: calculation for alternate allele frequency: p = (2*N_AA + N_AB) / 2N
# -> within genlight object (0,1,2): mean(dosage) / 2 is the alternate
# allele frequency

# function to calculate the alternate allele frequency per sampling site for one SNP
allele_freq_by_site <- function(genlight_obj, snp_id, run_table,
                                site_col = "geo_loc_name") {
  as.data.frame(as.matrix(genlight_obj)) |>
    rownames_to_column("Run") |>
    select(Run, dosage = all_of(snp_id)) |>
    left_join(run_table, by = "Run") |>
    group_by(across(all_of(site_col))) |>
    summarise(
      p         = mean(dosage, na.rm = TRUE) / 2,
      n         = sum(!is.na(dosage)),
      longitude = mean(longitude, na.rm = TRUE),
      latitude  = mean(latitude,  na.rm = TRUE),
      .groups = "drop"
    )
}

# retrieve one env. variable per sample site (as it is the same per site)
env_value_by_site <- function(recifs_sample, run_table, env_col,
                              site_col = "geo_loc_name") {
  recifs_sample |>
    rownames_to_column("Run") |>
    select(Run, env_value = all_of(env_col)) |>
    left_join(run_table, by = "Run") |>
    group_by(across(all_of(site_col))) |>
    summarise(
      env_value = dplyr::first(na.omit(env_value)),
      n         = dplyr::n(),                          # sample count at the site
      longitude = mean(longitude, na.rm = TRUE),
      latitude  = mean(latitude,  na.rm = TRUE),
      .groups = "drop"
    )
}

# per-sample dosage (0/1/2) for one SNP joined to that sample's own env. value
genotype_env_df <- function(genlight_obj, snp_id, recifs_sample, env_col) {
  dosage_df <- as.data.frame(as.matrix(genlight_obj)) |>
    rownames_to_column("Run") |>
    select(Run, dosage = all_of(snp_id))
  env_df <- recifs_sample |>
    rownames_to_column("Run") |>
    select(Run, env_value = all_of(env_col))
  dosage_df |>
    left_join(env_df, by = "Run") |>
    filter(!is.na(dosage), !is.na(env_value))
}

## little reminder on genotype coding:
# 0 = homozygous for the reference (or major) allele
# 1 = heterozygote (one reference/major allele, one alternative/minor allele)
# 2 = homozygous for the alternative (or minor) allele
## in relation to alternate allele frequency:
# -> more alternative allele = closer to 1, less alternate alleles = closer to 0

#==============================================================================
## 3) thesis summary figure: species title | map | GT-env (LFMM) | GT-env (CA)
#==============================================================================

#------------------------------------------------------------------------------
# 3.1) prepare data for plotting
#------------------------------------------------------------------------------

# linear rescale (e.g. sample counts -> circle radius range)
rescale_lin <- function(x, from, to) {
  if (diff(from) == 0) return(rep(mean(to), length(x)))
  to[1] + (x - from[1]) * (diff(to) / diff(from))
}

# build a circle polygon in plain lon/lat degrees (not geographically exact ->
# fine for a fixed-size marker glyph)
make_circle <- function(lon, lat, radius, id, n_points = 60) {
  angles <- seq(0, 360, length.out = n_points) * pi / 180
  data.frame(lon = lon + radius * sin(angles), lat = lat + radius * cos(angles), id = id)
}

# define env. variable and retreive target SNPs for that env. variable
selected_env <- "s_allDB_FE_sd" # check best_genes_lfmm (or ca)$env_var
selected_env_label <- "Iron Standard Deviation [mmol/m³]"    # full descriptive name -> spelled out once, in each figure's subtitle
selected_env_short_label <- "Iron SD [mmol/m³]"              
lfmm_locus <- best_genes_lfmm$locus[best_genes_lfmm$env_var == selected_env]

ca_locus   <- best_genes_ca$locus[best_genes_ca$env_var == selected_env]
target_snps_lfmm_sel <- get_target_snps(
  gene_id = lfmm_locus, env_var = selected_env, method = "LFMM",
  all_data = all_data, gene_results_field = "gene_lfmm", raw_results_field = "lfmm_results"
)
target_snps_ca_sel <- get_target_snps(
  gene_id = ca_locus, env_var = selected_env, method = "CA",
  all_data = all_data, gene_results_field = "gene_ca", raw_results_field = "ca_results"
)

#------------------------------------------------------------------------------
# 3.2) map panel: background env. circle + allele-freq circle
#------------------------------------------------------------------------------

# function to create summary plot
plot_summary_map <- function(env_df, af_df, chosen_method, same_snp,
                             world, xlim, ylim, env_label,
                             env_colours = c("#03045E", "#0077B6", "#00B4D8", "#90E0EF", "#CAF0F8"),
                             af_colours = c("white", "black"),
                             radius_frac = 0.03,              # smaller background env. circles
                             af_radius_base_frac = 0.045,     # allele-freq marker size, independent of the smaller background circle
                             af_radius_frac = c(0.16, 0.4),
                             n_legend_breaks = 3,             # circle size (N samples)
                             dark_threshold = 0.4,            # CHANGED: env-circle luminance below this value -> white af-circle outline, else black
                             land_fill = "grey55", land_colour = "grey40") {
  # env. var legend
  env_bar <- guide_colorbar(direction = "horizontal", title.position = "top",
                            barwidth = unit(4.5, "cm"), barheight = unit(0.35, "cm"),
                            theme = theme(legend.text = element_text(angle = 90, hjust = 1, vjust = 0.5)),
                            order = 1)
  # allele-frequency legend
  af_bar <- guide_colorbar(direction = "horizontal", title.position = "top",
                           barwidth = unit(4.5, "cm"), barheight = unit(0.35, "cm"),
                           order = 2)
  # background circle radius as a fraction of the map's longitude span -> keeps
  env_radius <- radius_frac * diff(xlim)
  af_radius_base <- af_radius_base_frac * diff(xlim)   
  site_df <- dplyr::inner_join(
    env_df |> dplyr::select(geo_loc_name, env_value, longitude, latitude),
    af_df  |> dplyr::select(geo_loc_name, p, n),
    by = "geo_loc_name"
  )
  n_range <- range(site_df$n)
  af_radius_for_n <- function(n) rescale_lin(sqrt(n), from = sqrt(n_range), to = af_radius_frac * af_radius_base)
  # background: one full circle per site, coloured by that site's env. value
  bg_polys <- purrr::map_dfr(seq_len(nrow(site_df)), function(i) {
    row <- site_df[i, ]
    make_circle(row$longitude, row$latitude, env_radius, id = paste0("env_", i))
  })
  bg_fill <- data.frame(id = paste0("env_", seq_len(nrow(site_df))), value = site_df$env_value)
  bg_polys <- dplyr::left_join(bg_polys, bg_fill, by = "id")
  # define wheter AF circle outline black or white (depend on env. colour/variation)
  env_ramp  <- grDevices::colorRamp(env_colours, space = "Lab")
  env_range <- range(site_df$env_value, na.rm = TRUE)
  env_norm  <- rescale_lin(site_df$env_value, from = env_range, to = c(0, 1))
  env_rgb   <- env_ramp(env_norm) / 255
  env_lum   <- 0.2126 * env_rgb[, 1] + 0.7152 * env_rgb[, 2] + 0.0722 * env_rgb[, 3]   # perceptual luminance
  outline_df <- data.frame(id = paste0("af_", seq_len(nrow(site_df))),
                           outline = ifelse(env_lum < dark_threshold, "white", "black"))
  # create plot for env. circles
  p <- ggplot() +
    geom_sf(data = world, fill = land_fill, colour = land_colour, linewidth = 0.3) +
    geom_polygon(data = bg_polys, aes(x = lon, y = lat, group = id, fill = value)) +
    scale_fill_gradientn(colours = env_colours, name = env_label, guide = env_bar)   
  # frequency of whichever method's SNP had the lower p-value
  af_polys <- purrr::map_dfr(seq_len(nrow(site_df)), function(i) {
    row <- site_df[i, ]
    make_circle(row$longitude, row$latitude, af_radius_for_n(row$n), id = paste0("af_", i))
  })
  af_fill <- data.frame(id = paste0("af_", seq_len(nrow(site_df))), value = site_df$p)
  af_polys <- dplyr::left_join(af_polys, af_fill, by = "id") |>
    dplyr::left_join(outline_df, by = "id")
  # legend title
  af_legend_name <- if (same_snp) "Allele Frequency" else paste0("Allele Frequency (", chosen_method, ")")
  p <- p + ggnewscale::new_scale_fill() +
    geom_polygon(data = af_polys, aes(x = lon, y = lat, group = id, fill = value, colour = outline),
                 linewidth = 0.8) +
    scale_fill_gradientn(colours = af_colours, name = af_legend_name,
                         limits = range(site_df$p, na.rm = TRUE),
                         labels = function(x) sprintf("%.1f", x), guide = af_bar) +
    scale_colour_identity() 
  # invisible layer to draw a single "N samples" size legend (empty circle, black outline, no fill)
  n_breaks <- round(seq(n_range[1], n_range[2], length.out = n_legend_breaks))
  p <- p +
    geom_point(data = site_df, aes(x = longitude, y = latitude, size = n), alpha = 0) +
    scale_size_continuous(name = "N samples", range = c(2, 6), breaks = n_breaks,
                          guide = guide_legend(override.aes = list(alpha = 1, shape = 21,
                                                                   fill = "white", colour = "black"),
                                               direction = "horizontal", title.position = "top",
                                               order = 3))
  # add annotations, north arrow, scale bar etc.
  p +
    coord_sf(xlim = xlim, ylim = ylim, expand = FALSE) +
    annotation_scale(location = "bl", width_hint = 0.2, style = "bar",
                     bar_cols = c("black", "white"), height = unit(0.4, "cm"),
                     text_cex = 1.1,                                        # CHANGED: only the label text is bigger (was the ggspatial default, ~0.7) -> bar itself back to its original size
                     unit_category = "metric", dist_unit = "km") +
    annotation_north_arrow(location = "tr", which_north = "true",
                           style = north_arrow_fancy_orienteering()) +
    labs(x = "Longitude", y = "Latitude",
         title = "Per-Site Allele Frequency & Environmental Variation",   # CHANGED: shortened ("and" -> "&")
         subtitle = if (same_snp) "Same SNP identified by LFMM and CA"
         else paste0("Different SNPs identified - ", chosen_method, " SNP plotted")) +   # CHANGED: shortened -> the reasoning (lower p-value) is explained in the methods/figure legend instead
    theme_minimal(base_size = 16) +   # CHANGED: bumped up slightly from 15 -> fewer rows per figure now, so more room per panel
    theme(panel.grid.major = element_line(colour = "grey90"),
          plot.title = element_text(face = "bold", size = 18, hjust = 0.5),   # CHANGED: bigger + centred
          plot.subtitle = element_text(size = 15, hjust = 0.5),               # CHANGED: bumped up (still below the 18pt title) + centred
          axis.title = element_text(size = 18),
          axis.text = element_text(size = 14),
          legend.title = element_text(size = 18),
          legend.text = element_text(size = 14),
          legend.position = "bottom", legend.box = "horizontal")
}

#------------------------------------------------------------------------------
# 3.3) boxplot panel: one method's genotype groups (3 boxes: 0/1/2)
#------------------------------------------------------------------------------

# function to create box plot for a SINGLE method
plot_summary_boxplot_single <- function(df, method_label, env_label, p_value,
                                        env_colours = c("#03045E", "#0077B6", "#00B4D8", "#90E0EF", "#CAF0F8"),
                                        geno_colours = c("0" = "grey90", "1" = "grey40", "2" = "grey10")) {
  y_rng <- range(df$env_value)
  # colour gradient background -> smooth raster gradient
  grad_colours <- grDevices::colorRampPalette(env_colours)(256)   # low -> high value ramp
  grad_matrix  <- matrix(rev(grad_colours), ncol = 1)             # row 1 = highest value, drawn at the top (ymax)
  # create plot
  ggplot() +
    annotation_raster(grad_matrix, xmin = -Inf, xmax = Inf, ymin = y_rng[1], ymax = y_rng[2],
                      interpolate = TRUE) +
    geom_boxplot(data = df, aes(x = dosage, y = env_value, group = dosage),
                 width = 0.5, outlier.shape = NA, fill = "white", alpha = 0.5, colour = "black") +
    # fixed seed -> when LFMM and CA share the same SNP
    geom_jitter(data = df, aes(x = dosage, y = env_value, fill = factor(dosage)),
                position = position_jitter(width = 0.12, height = 0, seed = 42),
                shape = 21, colour = "black", size = 2.2, stroke = 0.4) +
    # genotype colour legend -> guide order fixed
    scale_fill_manual(values = geno_colours, name = "Genotype", guide = guide_legend(order = 4)) +
    scale_x_continuous(breaks = c(0, 1, 2), labels = c("0", "1", "2")) +
    # strech env. gradient along x-axis (continuously)
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
    labs(title = paste0(method_label, ": Variation Across Genotypes"),                        
         subtitle = sprintf("p-value = %.4f", p_value),                                       
         x = "Genotype", y = env_label) +                                                     
    theme_classic(base_size = 16) +      # CHANGED: bumped up slightly from 15 -> fewer rows per figure now, so more room per panel
    theme(plot.title = element_text(face = "bold", size = 18),               
          plot.subtitle = element_text(size = 17),                          
          axis.title = element_text(size = 18),
          axis.text = element_text(size = 14),
          legend.title = element_text(size = 18),
          legend.text = element_text(size = 14),
          legend.position = c(0.06, 0.94),
          legend.justification = c(0, 1),
          legend.background = element_rect(fill = grDevices::adjustcolor("white", alpha.f = 0.7), colour = NA))
}

#------------------------------------------------------------------------------
# 3.4) combine rows to create the two summary figures
#------------------------------------------------------------------------------

# standardise every species' map to the SAME width:height ratio -> take the
# widest natural bounding box's ratio -> make others match!
target_map_ratio <- max(purrr::map_dbl(all_data, function(d) diff(d$xlim) / diff(d$ylim)))
pad_to_ratio <- function(xlim, ylim, target_ratio) {
  w <- diff(xlim); h <- diff(ylim)
  cur_ratio <- w / h
  if (cur_ratio < target_ratio) {
    new_w <- h * target_ratio
    pad <- (new_w - w) / 2
    xlim <- xlim + c(-pad, pad)
  } else if (cur_ratio > target_ratio) {
    new_h <- w / target_ratio
    pad <- (new_h - h) / 2
    ylim <- ylim + c(-pad, pad)
  }
  list(xlim = xlim, ylim = ylim)
}

# species titles: species name as small panel stacked above
make_title_strip <- function(label) {
  ggplot() +
    annotate("text", x = 0.5, y = 0, label = label, fontface = "bold.italic",
             size = 8, hjust = 0.5, vjust = 0) +
    xlim(0, 1) +
    coord_cartesian(clip = "off") +   # stop the large title text being clipped by the tiny title panel
    theme_void() +
    theme(plot.margin = margin(t = 10, r = 0, b = 5, l = 0))
}
# create rows
build_summary_row <- function(d, sp) {
  lfmm_snp <- target_snps_lfmm_sel$snp[target_snps_lfmm_sel$species == d$species]
  ca_snp   <- target_snps_ca_sel$snp[target_snps_ca_sel$species == d$species]
  p_lfmm   <- target_snps_lfmm_sel$p_value[target_snps_lfmm_sel$species == d$species] 
  p_ca     <- target_snps_ca_sel$p_value[target_snps_ca_sel$species == d$species]
  if (length(lfmm_snp) == 0 || is.na(lfmm_snp) || length(ca_snp) == 0 || is.na(ca_snp)) {
    message("Skipping '", sp, "' in the summary figure -> missing a target SNP for LFMM and/or CA")
    return(NULL)
  }
  same_snp <- identical(lfmm_snp, ca_snp)
  # map shows only ONE allele-frequency circle per site (method with SNP with lower p-value)
  chosen_method <- if (p_lfmm <= p_ca) "LFMM" else "CA"
  env_df     <- env_value_by_site(d$recifs_sample, d$run_table, selected_env)
  af_lfmm_df <- allele_freq_by_site(d$genlight, lfmm_snp, d$run_table)
  af_ca_df   <- allele_freq_by_site(d$genlight, ca_snp,   d$run_table)
  af_map_df  <- if (chosen_method == "LFMM") af_lfmm_df else af_ca_df 
  bounds     <- pad_to_ratio(d$xlim, d$ylim, target_map_ratio)
  map_panel <- plot_summary_map(env_df, af_map_df, chosen_method, same_snp,
                                world, xlim = bounds$xlim, ylim = bounds$ylim,
                                env_label = selected_env_short_label) 
  df_lfmm_gt <- genotype_env_df(d$genlight, lfmm_snp, d$recifs_sample, selected_env)
  df_ca_gt   <- genotype_env_df(d$genlight, ca_snp,   d$recifs_sample, selected_env)
  # two separate GT-env panels (one per method), each with its own rounded p-value
  box_panel_lfmm <- plot_summary_boxplot_single(df_lfmm_gt, method_label = "LFMM", env_label = selected_env_short_label, p_value = p_lfmm)  # CHANGED: short label + p-value
  box_panel_ca   <- plot_summary_boxplot_single(df_ca_gt,   method_label = "CA",   env_label = selected_env_short_label, p_value = p_ca)      # CHANGED: short label + p-value
  # LFMM and CA boxplots now stacked on top of each other
  row_panels <- (map_panel | (box_panel_lfmm / box_panel_ca)) +
    patchwork::plot_layout(guides = "collect", axes = "collect", axis_titles = "collect",
                           widths = c(0.6, 0.4)) &
    theme(legend.position = "bottom", legend.box = "horizontal", legend.justification = "center")
  make_title_strip(d$species) / row_panels + patchwork::plot_layout(heights = c(0.06, 1))
}
# combine all rows
summary_rows <- purrr::imap(all_data, build_summary_row)

# alphabetical order by species binomial: Amphiprion bicinctus, Epinephelus
# striatus, Pterois volitans, Siphamia tubifer
alpha_order  <- c("amphibic", "epistri", "ptevol", "siphtub")
summary_rows <- summary_rows[alpha_order] # reorder

# create figer IDs for split figure
figure_1_ids <- alpha_order[1:2]
figure_2_ids <- alpha_order[3:4]

# wrap a subset of rows into one annotated figure
build_summary_figure <- function(row_ids, panel_label = NULL) {
  rows <- purrr::compact(summary_rows[row_ids])
  title_suffix <- if (is.null(panel_label)) "" else paste0(" (", panel_label, ")")
  patchwork::wrap_plots(rows, ncol = 1) +
    patchwork::plot_annotation(
      title    = paste0("Repeated Adaptive Signal in Candidate Gene ", lfmm_locus,
                        " and Iron Concentration Variability", title_suffix),
      subtitle = paste0("Environmental variable: ", selected_env_label),
      theme = theme(plot.title = element_text(face = "bold", size = 26, hjust = 0.5, margin = margin(b = 4)),    # CHANGED: bumped up alongside the other titles
                    plot.subtitle = element_text(size = 17, hjust = 0.5, margin = margin(b = 15)))              # CHANGED: bumped up alongside the other subtitles
    )
}

# figure 1: Amphiprion bicinctus, Epinephelus striatus
summary_figure_1 <- build_summary_figure(figure_1_ids, "1/2")
# figure 2: Pterois volitans, Siphamia tubifer
summary_figure_2 <- build_summary_figure(figure_2_ids, "2/2")

# splitted summary figure (2 species per figure)
summary_figure_1
summary_figure_2
# format per figure: 1800 x 1800 (2 rows)

# all species summary figure (not for publications)
summary_figure_all <- build_summary_figure(alpha_order)
summary_figure_all
# format: 2000 x 4000 (4 rows)
