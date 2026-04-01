## SPDX-FileCopyrightText: 2026 Alkahest, Inc. <chronospd_data_access@alkahest.com>
##
## SPDX-License-Identifier: MIT

#' Export Results to Excel
#'
#' @param results joint_model_results object
#' @param file Path to output file
#' @export
export_results <- function(results, file = "joint_model_results.xlsx") {
  sheets <- list(
    summary = results$summaries |>
      dplyr::arrange(term, p_value),
    convergence = results$convergence,
    top_hits = results$summaries |>
      dplyr::filter(p_value < 0.05) |>
      dplyr::arrange(p_value) |>
      dplyr::mutate(
        hazard_ratio = exp(estimate),
        hr_lower = exp(lower_95),
        hr_upper = exp(upper_95),
        hr_95ci = glue::glue("{round(hazard_ratio, 2)} ({round(hr_lower, 2)}-{round(hr_upper, 2)})")
      ),
    failed_models = tibble::tibble(
      protein = results$failed_proteins,
      error = results$errors[results$failed_proteins]
    ),
    model_info = tibble::tibble(
      metric = names(results$data_info),
      value = unlist(results$data_info)
    )
  )

  writexl::write_xlsx(sheets, file)
  cli::cli_alert_success("Results exported to {file}")
}
