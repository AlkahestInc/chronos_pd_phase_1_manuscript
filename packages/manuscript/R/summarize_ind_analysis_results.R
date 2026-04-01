## SPDX-FileCopyrightText: 2026 Alkahest, Inc. <chronospd_data_access@alkahest.com>
##
## SPDX-License-Identifier: MIT

#' Create summary statistics table for protein and ratio analysis results
#'
#' @param results_table Tibble output from analyze_proteins_and_ratios function.
#' @param sig_threshold Numeric. P-value threshold for significance (default: 0.05).
#'
#' @return A list containing:
#'   - overall_summary: Overall summary statistics
#'   - ratio_summary: Summary for ratio analyses only
#'   - individual_summary: Summary for individual protein analyses
#'   - significant_ratios: Subset of significant ratio results
#'   - significant_individuals: Subset of significant individual protein results
#'   - failed_models: All analyses where models failed or didn't converge
#'
#' @examples
#' \dontrun{
#' summary <- summarize_analysis_results(results_table)
#' }
summarize_ind_analysis_results <- function(results_table, sig_threshold = 0.05) {
  individual_results <- results_table

  # Overall summary
  overall_summary <- tibble::tibble(
    total_analyses = nrow(results_table),
    total_individuals = nrow(individual_results),
    successful_models = sum(!is.na(results_table$estimate)),
    converged_models = sum(results_table$converged, na.rm = TRUE),
    significant_raw = sum(results_table$p_value < sig_threshold, na.rm = TRUE)
  )

  # Individual protein summary
  individual_summary <- tibble::tibble(
    total_proteins = nrow(individual_results),
    successful_models = sum(!is.na(individual_results$estimate)),
    converged_models = sum(individual_results$converged, na.rm = TRUE),
    significant_raw = sum(individual_results$p_value < sig_threshold, na.rm = TRUE),
    median_or = stats::median(individual_results$or, na.rm = TRUE),
    mean_n_obs = mean(individual_results$n_obs, na.rm = TRUE)
  )

  # Significant results
  significant_individuals <- individual_results |>
    dplyr::filter(p_value < sig_threshold) |>
    dplyr::select(
      analyte, protein_name, gene_symbol, estimate, or,
      or_lower_ci, or_upper_ci, p_value
    )

  # Failed models
  failed_models <- results_table |>
    dplyr::filter(is.na(estimate) | !converged) |>
    dplyr::select(analyte, protein_name, gene_symbol, converged)

  return(list(
    overall_summary = overall_summary,
    individual_summary = individual_summary,
    significant_individuals = significant_individuals,
    failed_models = failed_models
  ))
}
