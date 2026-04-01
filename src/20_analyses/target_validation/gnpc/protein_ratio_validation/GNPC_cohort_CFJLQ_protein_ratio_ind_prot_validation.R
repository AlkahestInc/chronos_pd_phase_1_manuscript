## SPDX-FileCopyrightText: 2026 Alkahest, Inc. <chronospd_data_access@alkahest.com>
##
## SPDX-License-Identifier: MIT

here::i_am("src/20_analyses/target_validation/gnpc/protein_ratio_validation/GNPC_cohort_CFJLQ_protein_ratio_ind_prot_validation.R")

input_dir <- "/home/rstudio"

# load GNPC annotation tibble
annotations <- glue::glue("{input_dir}/data_prep/annotation/annotation_tbl_v6.rds") |>
  readr::read_rds()

# load protein biomarkers
prot_ratio_tibble <- glue::glue(
  "{get_config()$input_dir}/name={name}/version={version}/dataset.parquet",
  name = "Chronos_phase1_protein_ratio_list_1",
  version = 29452
) |>
  arrow::open_dataset() |>
  dplyr::collect()

# probe names from somalogic starts with `seq`, from olink starts with `OID`, and alamarbio directly given gene symbols
prot_ratio_tibble_list <- list(
  chronos_phase1_protein_ratio_list_1 = prot_ratio_tibble
)

# check/create output path
output_dir <- glue::glue("{get_config()$checkpoint_dir}/GNPC_validation/model/protein_ratio_validation/ind_protein_validation_result")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Build "other" analyte mapping only if missing
# other_analyte_mapping Data frame for other platform mappings.
mapping_tibbles <- prot_ratio_tibble_list |>
  build_mapping_tibbles(
    protein1_col = "protein1",
    protein2_col = "protein2",
    uniprot1_col = "uniprot1",
    uniprot2_col = "uniprot2",
    gene_symbol1_col = "gene_symbol1",
    gene_symbol2_col = "gene_symbol2"
  )

other_analyte_mapping <- mapping_tibbles |>
  purrr::reduce(dplyr::bind_rows) |>
  dplyr::distinct(.keep_all = TRUE)

for (selected_contributor_code in c("C", "F", "J_plasma", "L", "Q_plasma")) {
  cat("For GNPC cohort : ", selected_contributor_code, "\n")

  if (selected_contributor_code == "J_plasma") {
    matching_ratio_str <- "1to1"
  } else {
    matching_ratio_str <- "2to1"
  }
  if (selected_contributor_code %in% c("C", "L", "Q_plasma")) {
    visit_str <- ""
  } else {
    visit_str <- "_last_visit" # "_first_visit"
  }

  # load gnpc protein intensity and other data
  metadata <- glue::glue("{get_config()$checkpoint_dir}/data_prep/gnpc_cohort_cleanup_matched_data/cohort_{selected_contributor_code}{visit_str}_outliers_rm_{matching_ratio_str}_matched_grp_filtered_metadata.rds") |>
    readr::read_rds()
  protein_intensity_log2 <- glue::glue("{get_config()$checkpoint_dir}/data_prep/gnpc_cohort_cleanup_matched_data/cohort_{selected_contributor_code}{visit_str}_outliers_rm_{matching_ratio_str}_matched_grp_filtered_protein_intensity_log2.rds") |>
    readr::read_rds()

  results <- tibble::tibble(
    analyte = character(),
    protein_name = character(),
    gene_symbol = character(),
    estimate = numeric(),
    or_lower_ci = numeric(),
    or_upper_ci = numeric(),
    p_value = numeric(),
    n_obs = integer(),
    converged = logical()
  )

  proteins <- unique(c(prot_ratio_tibble$protein1, prot_ratio_tibble$protein2))
  analytes <- c()
  gene_symbols <- c()
  for (i in seq_along(proteins)) {
    proc_analytes <- unique(process_protein_ratio_prob_str(proteins[i],
      mapping_cross_olink_platforms = FALSE,
      mapping_cross_somalogic_platforms = FALSE
    ))
    if (length(proc_analytes) == 0 && is.character(proc_analytes)) {
      proc_analytes <- c(proteins[i])
    }

    proc_gene_symbols <- unique(get_gene_symbol_by_prob(proteins[i],
      mapping_cross_olink_platforms = FALSE,
      mapping_cross_somalogic_platforms = FALSE
    ))
    if (length(proc_analytes) == 0 && is.character(proc_analytes)) {
      proc_gene_symbols <- c(proteins[i])
    }
    analytes <- c(analytes, proc_analytes)
    gene_symbols <- c(gene_symbols, proc_gene_symbols)

    analyte_gene_symbols <- expand.grid(
      analyte = proc_analytes,
      gene_symbol = proc_gene_symbols,
      stringsAsFactors = FALSE
    )
    new_res <- purrr::pmap_dfr(
      analyte_gene_symbols,
      ~ analyze_single_protein(proteins[i],
        ..1,
        ..2,
        analysis_data = metadata,
        analysis_prot = protein_intensity_log2
      )
    )
    results <- results |>
      tibble::add_row(!!!new_res)
  }

  # Add additional calculated columns
  results <- results |>
    dplyr::mutate(or = exp(estimate)) |>
    dplyr::arrange(p_value)

  summaries <- summarize_ind_analysis_results(results_table = results)

  # filter out unmatched or unconverged rows (NA rows)
  results <- results |>
    dplyr::filter(converged)

  output_file_path <- glue::glue("{output_dir}/cohort_{selected_contributor_code}{visit_str}_outliers_rm_{matching_ratio_str}_matched_grp_3tbl_protein_ratio_ind_protein_validation_results.csv")
  dir.create(dirname(output_file_path), showWarnings = FALSE, recursive = TRUE)
  results |>
    readr::write_csv(output_file_path)
  cat("Output individual protein validation result table for GNPC cohort : ", selected_contributor_code, "; results dims (", dim(results), ")\n")

  output_summary_significant_path <- glue::glue("{output_dir}/cohort_{selected_contributor_code}{visit_str}_outliers_rm_{matching_ratio_str}_matched_grp_3tbl_protein_ratio_ind_protein_validation_result_summary_significant_ind.csv")
  dir.create(dirname(output_summary_significant_path), showWarnings = FALSE, recursive = TRUE)
  summaries$significant_individuals |>
    readr::write_csv(output_summary_significant_path)
  cat("Output significant individual protein table for GNPC cohort : ", selected_contributor_code, "; results dims (", dim(summaries$significant_individuals), ")\n")

  output_summary_stats_path <- glue::glue("{output_dir}/cohort_{selected_contributor_code}{visit_str}_outliers_rm_{matching_ratio_str}_matched_grp_3tbl_protein_ratio_ind_protein_validation_summary_stats.csv")
  dir.create(dirname(output_summary_stats_path), showWarnings = FALSE, recursive = TRUE)
  summaries$individual_summary |>
    readr::write_csv(output_summary_stats_path)

  print(summaries$individual_summary)
  cat("Output individual protein validation summary stats table for GNPC cohort : ", selected_contributor_code, "\n")
  cat("Done for GNPC cohort : ", selected_contributor_code, "\n")
}
