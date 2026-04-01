## SPDX-FileCopyrightText: 2026 Alkahest, Inc. <chronospd_data_access@alkahest.com>
##
## SPDX-License-Identifier: MIT

#' Create joint_model_results object
#' @keywords internal
joint_model_results <- function(results, proteins, spec, data) {
  # Extract successful models
  successful <- purrr::keep(results, ~ .x$status == "success")
  failed <- purrr::keep(results, ~ .x$status == "failed")

  structure(
    list(
      models = stats::setNames(
        purrr::map(successful, "model"),
        purrr::map_chr(successful, "protein")
      ),
      summaries = purrr::map_dfr(successful, "summary"),
      convergence = purrr::map_dfr(successful, ~ tibble::tibble(
        protein = .x$protein,
        converged = .x$convergence$converged,
        max_rhat = .x$convergence$max_rhat,
        min_ess = .x$convergence$min_ess
      )),
      failed_proteins = purrr::map_chr(failed, "protein"),
      errors = stats::setNames(
        purrr::map_chr(failed, "error"),
        purrr::map_chr(failed, "protein")
      ),
      spec = spec,
      data_info = list(
        n_subjects = data$n_subjects,
        n_observations = data$n_observations,
        n_proteins = length(proteins),
        n_successful = length(successful),
        n_failed = length(failed)
      ),
      call = match.call()
    ),
    class = c("joint_model_results", "list")
  )
}
