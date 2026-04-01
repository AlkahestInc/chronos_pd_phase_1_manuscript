## SPDX-FileCopyrightText: 2026 Alkahest, Inc. <chronospd_data_access@alkahest.com>
##
## SPDX-License-Identifier: MIT

#' Print method for joint_model_results
#' @export
print.joint_model_results <- function(x, ...) {
  cli::cli_h1("Joint Model Results")

  info <- x$data_info
  cli::cli_alert_success("{info$n_successful}/{info$n_proteins} models fitted successfully")

  if (info$n_failed > 0) {
    cli::cli_alert_warning("{info$n_failed} models failed")
  }

  # Convergence summary
  conv <- x$convergence
  n_converged <- sum(conv$converged, na.rm = TRUE)
  cli::cli_alert_info("{n_converged}/{info$n_successful} models converged (Rhat < 1.05)")

  # Top associations
  if (nrow(x$summaries) > 0) {
    top <- x$summaries |>
      dplyr::filter(term == "value") |>
      dplyr::arrange(p_value) |>
      dplyr::slice_head(n = 5)

    if (nrow(top) > 0) {
      cli::cli_h2("Top 5 Associations")
      cli::cli_ul(glue::glue("{top$protein}: HR = {round(exp(top$estimate), 2)} (p = {format.pval(top$p_value, digits = 2)})"))
    }
  }

  invisible(x)
}
