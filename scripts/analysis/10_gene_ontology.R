#==============================================================================
### tenth analysis - gene ontology analysis for picmin outputs
#==============================================================================

# This script runs a gene set enrichment analysis using SetRank on the PicMin
# outputs, for every environmental variable that has at least one gene passing
# the pooled q-value threshold (q_threshold, currently 0.25). Runs all three
# GO domains (BP, MF, CC); only BP is used in the summary figure (script 12),
# MF and CC results are kept for later use.
# Inputs:
# 1) picmin output, results list (CA and LFMM)
# 2) danio rerio uniprot reference file: Danio_rerio.GRCz11.116.uniprot.tsv
# 3) ENSEMBL <-> UniProt linker: uniprotkb_AND_model_organism_7955_2026_04_29.tsv
# Output:
# 1. gene ontology results
# 2. GO term summary table
# 3. two GO enrichment heatmaps, one per method (LFMM, CA)

# set correct working directory:
setwd("C:/Users/Admin/Documents/Uni_Zuerich/Masters_Thesis/CoralReef_fish_local")

# Author: Oliver M. Selmoni, adjusted by Naroa Schweizer
# last update: 20.08.26

# load libraries
library(purrr)           # iterating over lists/vectors without writing for loops
library(dplyr)           # data wrangling (filter, summarize, join)
library(SetRank)         # gene set enrichment method
library(ggplot2)         # the two GO heatmaps (LFMM, CA) built in section 3

# NOTE: SetRank is a gene set enrichment method. Given a ranked list of
# candidate genes, it tests which GO terms (BP/MF/CC) or pathways show up
# more often among them than expected by chance, correcting for the fact
# that gene sets overlap heavily. This is a different q-value to PicMin's
# pooled_q: pooled_q scores individual genes, SetRank's output scores GO
# terms, so one cannot be derived from the other.

#------------------------------------------------------------------------------
## define paths and load files
#------------------------------------------------------------------------------

# define path to UniProt GO annotations for Danio rerio
uniprot_path <- "intermed_outputs/danio_rerio/uniprotkb_AND_model_organism_7955_2026_04_29.tsv"
# load annotation file
UNIPROT <- read.table(uniprot_path, sep = "\t", quote = "", header = TRUE)
rownames(UNIPROT) <- UNIPROT$Entry.Name

# define path to ENSEMBL <-> UniProt linker
ens_up_path <- "intermed_outputs/danio_rerio/Danio_rerio.GRCz11.116.uniprot.tsv"
# load UniProt <-> ENSEMBL linker
ENS_UP_full <- read.table(ens_up_path, header = TRUE, sep = "\t", quote = "")

# load picmin outputs for ca and lfmm
load("intermed_outputs/picmin/picmin_results_lfmm.rda"); lfmm_picmin <- picmin_results
load("intermed_outputs/picmin/picmin_results_ca.rda");   ca_picmin   <- picmin_results

# create output directory
# dir.create('intermed_outputs/setrank', showWarnings = FALSE)
# define output path to save etc.
out_dir <- "intermed_outputs/setrank"

# define number of species in specific PicMin run (here 4 -> A. bicinctus, E. striatus, P. volitans, S. tubifer)
n_species_required <- 4

# define q-value threshold for picmin
q_threshold <- 0.25

#=============================================================================
## 1) prepare inputs to run ontology: picmin gene universe
#=============================================================================

#-----------------------------------------------------------------------------
## 1.1) prepare GO annotation table
#-----------------------------------------------------------------------------

# function to build one GO annotation table (one gene <-> one GO term per
# row); written once and reused for BP (biological process), MF (molecular
# function) and CC (cellular component)
build_annotation_table <- function(ens_up, go_col, db_name) {
  do.call(rbind, by(ens_up, ens_up$gene_stable_id, function(g) {
    
    # take all the annotations associated with the gene
    allAnn <- strsplit(g[[go_col]], "; ")
    
    # keep only unique annotations
    allAnn <- unique(unlist(allAnn))
    
    # split GO term ID from GO term name
    termID   <- gsub("(.*) \\[(GO.*)\\]", "\\2", allAnn)
    termName <- gsub("(.*) \\[(GO.*)\\]", "\\1", allAnn)
    
    # prepare other columns for output
    geneID <- unique(g$gene_stable_id)  # gene ID
    description <- "-"                  # we leave this empty
    dbName <- db_name                   # name of GO database
    
    # set to NA if fields are missing
    if (length(termID) == 0)   termID   <- NA
    if (length(termName) == 0) termName <- NA
    
    # return as dataframe
    data.frame(geneID, termID, termName, description, dbName)
  }))
}

#-----------------------------------------------------------------------------
## 1.2) create picmin gene universe
#-----------------------------------------------------------------------------

# function to build & save the BP/MF/CC gene set collections for one method's
# PicMin gene universe (all loci PicMin ever tested for that method, across
# every env variable)
build_gene_universe <- function(picmin_path, method_label) {
  gsc_paths <- setNames(
    file.path(out_dir, paste0(c("BP", "MF", "CC"), "_gsc_", method_label, ".rda")),
    c("BP", "MF", "CC")
  )
  # gene <-> GO term annotation table, kept so enriched terms can be traced
  # back to the Ensembl gene(s) behind them later (section 3)
  annotation_path <- file.path(out_dir, paste0("annotation_", method_label, ".rda"))
  
  # skip rebuilding if everything already exists on disk from a previous run
  if (all(file.exists(c(gsc_paths, annotation_path)))) {
    message("BP/MF/CC gene set collections + annotation table for ", method_label, " already exist, skipping rebuild")
    return(invisible(NULL))
  }
  
  # load picmin results
  load(picmin_path)  # loads 'picmin_results'
  
  # keep only annotations in PicMin results
  ens_up <- ENS_UP_full[ENS_UP_full$gene_stable_id %in% picmin_results[[1]]$locus, ]
  
  # attach GO annotations to Ensembl identifiers
  ens_up$MF <- UNIPROT[ens_up$xref, "Gene.Ontology..molecular.function."]  # add molecular function ontologies
  ens_up$BP <- UNIPROT[ens_up$xref, "Gene.Ontology..biological.process."]  # add biological process ontologies
  ens_up$CC <- UNIPROT[ens_up$xref, "Gene.Ontology..cellular.component."]  # add cellular component ontologies
  
  # format the annotation table for SetRank (one gene <-> one GO annotation per line)
  annotationTABLEbp <- build_annotation_table(ens_up, "BP", "BP")
  annotationTABLEmf <- build_annotation_table(ens_up, "MF", "MF")
  annotationTABLEcc <- build_annotation_table(ens_up, "CC", "CC")
  
  # construct gene set collections
  BP_gsc <- buildSetCollection(annotationTABLEbp)
  MF_gsc <- buildSetCollection(annotationTABLEmf)
  CC_gsc <- buildSetCollection(annotationTABLEcc)
  
  # save the gene set collections
  save(BP_gsc, file = gsc_paths[["BP"]])
  save(MF_gsc, file = gsc_paths[["MF"]])
  save(CC_gsc, file = gsc_paths[["CC"]])
  
  # save the gene <-> GO term links (all three domains) for gene-level tracing
  annotation_all <- dplyr::bind_rows(annotationTABLEbp, annotationTABLEmf, annotationTABLEcc)
  save(annotation_all, file = annotation_path)
}

options(mc.cores = 6)

#-----------------------------------------------------------------------------
## 1.3) run gene universe creation for both methods
#-----------------------------------------------------------------------------

# build one gene universe per method (LFMM and CA)
build_gene_universe("intermed_outputs/picmin/picmin_results_lfmm.rda", "lfmm")
build_gene_universe("intermed_outputs/picmin/picmin_results_ca.rda",   "ca")

#=============================================================================
## 2) prepare inputs to run ontology: select environments, do setranks
#=============================================================================

#-----------------------------------------------------------------------------
## 2.1) select "best" environmental variables
#-----------------------------------------------------------------------------

# function to select environments with at least one gene below q_threshold
has_best_gene <- function(env_var, picmin_obj) {
  df <- picmin_obj[[env_var]]
  nrow(dplyr::filter(df, n_est == n_species_required, pooled_q < q_threshold)) > 0
}

#-----------------------------------------------------------------------------
## 2.2) create SetRank
#-----------------------------------------------------------------------------

# function to run SetRank (for BP/MF/CC) for one env. variable and export results
run_setrank_for_env <- function(env_var, picmin_obj, method_label, BP_gsc, MF_gsc, CC_gsc) {
  
  # select PicMin q-values, and use them to rank genes
  PM_RES <- picmin_obj[[env_var]]
  
  # keep only results below the pooled q-value threshold (q_threshold)
  PM_RES <- PM_RES[PM_RES$pooled_q < q_threshold, ]
  
  # rank genes by pooled q-value
  TOPgenes <- PM_RES$locus[order(PM_RES$pooled_q)]
  
  # run setrank for one GO domain and export the results
  run_one <- function(gsc, db_name) {
    out_path      <- file.path(out_dir, paste0("SetRankResults_", method_label, "_", env_var, "_", db_name))
    pathways_file <- paste0(out_path, "_pathways.txt")
    empty_marker  <- paste0(out_path, "_EMPTY")  # marks "ran, found nothing significant"
    
    # if this env. variable was already processed in a previous run, don't
    # run SetRank again, just read back the q-value that's already there
    if (file.exists(pathways_file)) {
      message("  ", db_name, " results for ", env_var, " already exist, reading cached result")
      return(read.table(pathways_file, header = TRUE, sep = "\t"))
    }
    if (file.exists(empty_marker)) {
      message("  ", db_name, " for ", env_var, " previously had 0 significant gene sets, skipping")
      return(NULL)
    }
    
    # run setrank
    net <- setRankAnalysis(TOPgenes, gsc, use.ranks = TRUE, setPCutoff = 0.1, fdrCutoff = 0.1)
    
    # if 0 gene sets come out significant, SetRank's exportSingleResult()
    # fails trying to order a column that doesn't exist; skip the export,
    # leave a marker behind, and return NULL instead
    if (igraph::vcount(net) == 0) {
      message("  no significant ", db_name, " gene sets for ", env_var)
      file.create(empty_marker)
      return(NULL)
    }
    
    # export network and read back the pathways table
    exportSingleResult(network = net, selectedGenes = TOPgenes, collection = gsc,
                       networkName = out_path, IDConverter = NULL)
    read.table(pathways_file, header = TRUE, sep = "\t")
  }
  
  # run setrank - BP, MF, CC
  list(
    BP = run_one(BP_gsc, "BP"),
    MF = run_one(MF_gsc, "MF"),
    CC = run_one(CC_gsc, "CC")
  )
}

# function to run SetRank over one full GEA method, across every env. variable
run_setrank_for_method <- function(picmin_path, gsc_label, method_label) {
  # load picmin results
  load(picmin_path)  # loads 'picmin_results'
  
  # load gene set collections
  load(file.path(out_dir, paste0("BP_gsc_", gsc_label, ".rda")))
  load(file.path(out_dir, paste0("MF_gsc_", gsc_label, ".rda")))
  load(file.path(out_dir, paste0("CC_gsc_", gsc_label, ".rda")))
  
  # select environments worth running (at least one gene below q_threshold)
  env_vars <- names(picmin_results)
  best_env_vars <- env_vars[sapply(env_vars, has_best_gene, picmin_obj = picmin_results)]
  message(method_label, ": running SetRank on ", length(best_env_vars), " env variable(s): ",
          paste(best_env_vars, collapse = ", "))
  
  # run setrank for every selected env. variable
  purrr::map(best_env_vars, run_setrank_for_env, picmin_obj = picmin_results,
             method_label = method_label, BP_gsc = BP_gsc, MF_gsc = MF_gsc, CC_gsc = CC_gsc) |>
    setNames(best_env_vars)
}

#-----------------------------------------------------------------------------
## 2.3) run SetRank for both methods
#-----------------------------------------------------------------------------

# run for both methods
results_lfmm <- run_setrank_for_method("intermed_outputs/picmin/picmin_results_lfmm.rda", "lfmm", "LFMM")
results_ca   <- run_setrank_for_method("intermed_outputs/picmin/picmin_results_ca.rda",   "ca",   "CA")

#=============================================================================
## 3) GO enrichment summary + heatmaps, one per method (LFMM, CA)
#=============================================================================

#-----------------------------------------------------------------------------
## 3.1) shared row order for both plots (env. variables with a hit in either method)
#-----------------------------------------------------------------------------

# n hits below q_threshold per env. variable, combined across methods -> descending order
hit_counts <- function(picmin_obj) sapply(picmin_obj, function(df) sum(df$pooled_q < q_threshold, na.rm = TRUE))
all_vars <- union(names(lfmm_picmin), names(ca_picmin))
combined_counts <- setNames(rep(0, length(all_vars)), all_vars)
for (counts in list(hit_counts(lfmm_picmin), hit_counts(ca_picmin))) combined_counts[names(counts)] <- combined_counts[names(counts)] + counts
env_order <- names(sort(combined_counts[combined_counts > 0], decreasing = TRUE))

#-----------------------------------------------------------------------------
## 3.2) run/load GO term summary, incl. the Ensembl genes behind each term
#-----------------------------------------------------------------------------

go_summary_path <- file.path(out_dir, "picmin_go_terms.csv")

if (file.exists(go_summary_path)) {
  
  message("GO term summary already exists, skipping rebuild: ", go_summary_path)
  go_summary <- read.csv(go_summary_path, stringsAsFactors = FALSE)
  
} else {
  
  # genes below q_threshold per env. variable, for ONE method
  hits_for_method <- function(picmin_obj) {
    dplyr::bind_rows(Map(function(df, env_name) {
      df <- df[df$pooled_q < q_threshold, , drop = FALSE]
      if (nrow(df) == 0) return(NULL)
      data.frame(env_var = env_name, locus = df$locus, stringsAsFactors = FALSE)
    }, picmin_obj, names(picmin_obj)))
  }
  
  # one env. variable / one ontology -> GO id + description + q-value + the
  # candidate genes
  run_go_enrichment <- function(env_name, hits_df, go_collection, annotation) {
    genes <- unique(hits_df$locus[hits_df$env_var == env_name])
    if (length(genes) == 0) return(NULL)
    
    net <- SetRank::setRankAnalysis(genes, go_collection, use.ranks = FALSE,
                                    setPCutoff = 0.01, fdrCutoff = 0.05)
    if (igraph::vcount(net) == 0) return(NULL)
    
    res <- igraph::as_data_frame(net, what = "vertices")
    genes_ensembl <- unname(sapply(res$name, function(id) {
      paste(sort(intersect(genes, annotation$geneID[annotation$termID == id])), collapse = "; ")
    }))
    
    data.frame(env_var = env_name, go_id = res$name, go_term = res$description,
               fdr = res$adjustedPValue, genes_ensembl = genes_ensembl, stringsAsFactors = FALSE)
  }
  
  options(mc.cores = 1)  # serial, matches the original setting for this call
  
  # run BP/CC/MF for one method, using that method's own gsc + annotation table
  build_method_summary <- function(picmin_obj, method_label, gsc_label) {
    load(file.path(out_dir, paste0("BP_gsc_", gsc_label, ".rda")))
    load(file.path(out_dir, paste0("MF_gsc_", gsc_label, ".rda")))
    load(file.path(out_dir, paste0("CC_gsc_", gsc_label, ".rda")))
    load(file.path(out_dir, paste0("annotation_", gsc_label, ".rda")))  # loads 'annotation_all'
    
    hits_df <- hits_for_method(picmin_obj)
    out <- purrr::imap_dfr(list(BP = BP_gsc, CC = CC_gsc, MF = MF_gsc), function(gsc, ontology_name) {
      res <- dplyr::bind_rows(lapply(env_order, run_go_enrichment,
                                     hits_df = hits_df, go_collection = gsc, annotation = annotation_all))
      if (nrow(res) > 0) res$ontology <- ontology_name
      res
    })
    out$method <- method_label
    out
  }
  
  go_summary <- dplyr::bind_rows(
    build_method_summary(lfmm_picmin, "LFMM", "lfmm"),
    build_method_summary(ca_picmin,   "CA",   "ca")
  )
  
  write.csv(go_summary, go_summary_path, row.names = FALSE)
  
}

#-----------------------------------------------------------------------------
## 3.3) two GO heatmaps, one per method
#-----------------------------------------------------------------------------

strip_env_prefix <- function(x) sub("^s_allDB_?", "", x)

row_grid_theme <- theme(
  panel.spacing.y    = unit(0, "lines"),
  panel.border       = element_rect(colour = "grey65", fill = NA, linewidth = 0.3),
  panel.grid.major.y = element_blank(),
  panel.grid.minor   = element_blank(),
  axis.text.y  = element_blank(),
  axis.ticks.y = element_blank()
)

# same look as the old pooled panel C plot, just filtered to one method
plot_go_heatmap <- function(method_label) {
  plot_data <- dplyr::filter(go_summary, method == method_label)
  
  # if this method has 0 significant terms anywhere, there's nothing to plot
  if (nrow(plot_data) == 0) {
    message("no significant GO terms for ", method_label, ", skipping plot")
    return(invisible(NULL))
  }
  
  plot_data$env_var  <- factor(plot_data$env_var, levels = env_order)
  plot_data$ontology <- factor(plot_data$ontology, levels = c("BP", "CC", "MF"))
  plot_data$go_term  <- factor(plot_data$go_term, levels = sort(unique(plot_data$go_term)))
  
  ggplot(plot_data, aes(x = go_term, y = "row", fill = fdr)) +
    geom_tile(colour = "grey85") +
    facet_grid(rows = vars(env_var), cols = vars(ontology), switch = "both", drop = FALSE,
               scales = "free_x",
               labeller = labeller(env_var = strip_env_prefix)) +
    scale_fill_gradientn(colours = c("#03045E", "#0077B6", "#00B4D8", "#90E0EF", "#CAF0F8"),
                         limits = c(0, 0.05), name = "GO q-value") +
    labs(x = "GO term, grouped by ontology (BP = biological process, CC = cellular component, MF = molecular function)",
         y = NULL,
         title = paste0("Gene Ontology Terms across Environmental Variables for PicMin Genes - ",
                        method_label, " (q-value < 0.25)")) +
    theme_minimal(base_size = 11) +
    theme(strip.placement = "outside",
          strip.text.y.left = element_text(angle = 0, hjust = 1),
          strip.text.x.bottom = element_text(face = "bold"),
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 9),
          panel.grid.major.x = element_line(colour = "grey85", linewidth = 0.3),
          panel.grid.major.y = element_blank(),
          panel.grid.minor = element_blank()) +
    row_grid_theme
}

plot_go_lfmm <- plot_go_heatmap("LFMM")
plot_go_ca   <- plot_go_heatmap("CA")

# display
plot_go_lfmm
# save: 800 x 500
plot_go_ca
# save: 800 x 500
