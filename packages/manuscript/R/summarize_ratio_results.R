## SPDX-FileCopyrightText: 2026 Alkahest, Inc. <chronospd_data_access@alkahest.com>
##
## SPDX-License-Identifier: MIT

#' Summarize protein ratio modeling results
#'
#' Computes high-level summary statistics, extracts significant ratios,
#' and identifies failed or non-converged models from a results table.
#'
#' @param results_table Data frame or tibble containing model results.
#' @param sig_threshold Numeric. Significance threshold for p-values. Default to 0.05
#'
#' @return A named list with elements:
#' \describe{
#'   \item{summary_stats}{One-row tibble of aggregate statistics}
#'   \item{significant_ratios}{Subset of significant model results}
#'   \item{failed_models}{Subset of failed or non-converged models}
#' }
summarize_ratio_results <- function(results_table, sig_threshold = 0.05) {
  summary_stats <- tibble::tibble(
    total_pairs = nrow(results_table),
    successful_models = sum(!is.na(results_table$estimate)),
    converged_models = sum(results_table$converged, na.rm = TRUE),
    significant_raw = sum(results_table$p_value < sig_threshold, na.rm = TRUE),
    median_or = stats::median(results_table$or, na.rm = TRUE),
    mean_n_obs = mean(results_table$n_obs, na.rm = TRUE)
  )
  significant_ratios <- results_table |>
    dplyr::filter(p_value < sig_threshold) |>
    dplyr::select(
      protein_pair, mapped_analyte_1, mapped_analyte_2, prot_1_uniprot,
      prot_2_uniprot, gene_symbol1, gene_symbol2, estimate, or,
      or_lower_ci, or_upper_ci, p_value
    )
  failed_models <- results_table |>
    dplyr::filter(is.na(estimate) | !converged) |>
    dplyr::select(
      protein_pair, mapped_analyte_1, mapped_analyte_2, prot_1_uniprot,
      prot_2_uniprot, gene_symbol1, gene_symbol2, converged
    )
  return(list(
    summary_stats = summary_stats,
    significant_ratios = significant_ratios,
    failed_models = failed_models
  ))
}
