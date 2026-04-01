## SPDX-FileCopyrightText: 2026 Alkahest, Inc. <chronospd_data_access@alkahest.com>
##
## SPDX-License-Identifier: MIT

here::i_am("src/20_analyses/sepa/SEPA_classic.R")

# ==============================================================================
# SEPA (Stepwise Enrichment Pathway Analysis) using modified PA_analysis
# ==============================================================================
# Author: T.Nanasi
# Date: 2025.11
# Description: Pathway enrichment analysis on Chronos proteomics data
#              (Somalogic, Olink, Alamar platforms) using GO, KEGG, and Reactome

# ==============================================================================
# ENVIRONMENT SETUP
# ==============================================================================
script_start_time <- Sys.time()

# Check R version for annotation package compatibility
if (as.numeric(version$major) < 4 | as.numeric(version$minor) < 4.1) {
  print("consider moving to R 4.4.1+ -- mappings in PA packages tend to get obsolete if not using the most recent R")
}

# ==============================================================================
# FILEPATHS - modify these as needed
# ==============================================================================
OUTPUT_DIR_PA <- here::here("src/20_analyses/sepa/PA")
OUTPUT_FILE_PA_READY <- here::here(glue::glue("{get_config()$checkpoint_dir}/chronos_PA_ready.rds"))
dir.create(dirname(OUTPUT_FILE_PA_READY), showWarnings = FALSE, recursive = TRUE)

# ==============================================================================
# DATA PREPARATION
# ==============================================================================

source(
  here::here("src/20_analyses/sepa/SEPA_data_preparation.R")
)

# ----------------------------------------------------------------------------
# Save prepared data for posterity
# ----------------------------------------------------------------------------

# Compatibility: restore original variable names for saved object
dats <- raw_datasets # compatibility
dats1 <- cleaned_datasets # compatibility
m <- annotation_maps # compatibility
paranks <- pathway_ranks # compatibility

X <- list()
X$dats <- dats
X$dats1 <- dats1
X$m <- m
X$paranks <- paranks

X |>
  readr::write_rds(OUTPUT_FILE_PA_READY)

rm(dats, dats1, m, paranks)

# ==============================================================================
# RUN PATHWAY ANALYSIS
# ==============================================================================

sepa_depths <- c(seq(10, 500, by = 10)) # SEPA depth levels

pracma::tic()

# Loop across condition/directionality constellations
# (i.e., preDx upregulated, all samples any direction, etc.)
for (q in seq_along(pathway_ranks)) {
  condition_name <- names(pathway_ranks)[[q]]
  print(paste("running PA for", condition_name))

  output_filename <- glue::glue("{get_config()$checkpoint_dir}/SEPA_UNIPROT_{condition_name}.Rdata")
  dir.create(dirname(output_filename), showWarnings = FALSE, recursive = TRUE)

  entity_rank_df <- pathway_ranks[[q]]$ranks
  background_ids <- pathway_ranks[[q]]$background

  # Create entity lists for each SEPA depth
  entity_sets <- list()
  for (k in 1:length(sepa_depths)) {
    entity_sets[[k]] <- entity_rank_df$entity[which(entity_rank_df$rank <= sepa_depths[k])]
  }
  entity_sets <- unique(entity_sets)

  # Run pathway enrichment analysis
  PA <- runPAs_parallel(
    entity_sets,
    background_ids,
    id_type = pathway_ranks[[q]]$universe
  )

  # Compatibility: restore original variable names for saved output
  used.ranks <- entity_rank_df # compatibility
  used.filters <- condition_name # compatibility
  used.topns <- sepa_depths # compatibility

  save(
    list = c("PA", "used.ranks", "used.filters", "used.topns"),
    file = output_filename
  )

  rm(entity_rank_df, background_ids, entity_sets, condition_name, output_filename)
  rm(used.ranks, used.filters, used.topns)

  pracma::toc()
}

print(paste("Total runtime:", round(difftime(Sys.time(), script_start_time, units = "mins"), 2), "minutes"))
