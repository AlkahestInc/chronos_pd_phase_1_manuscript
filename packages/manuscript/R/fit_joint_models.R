## SPDX-FileCopyrightText: 2026 Alkahest, Inc. <chronospd_data_access@alkahest.com>
##
## SPDX-License-Identifier: MIT

#' Fit Joint Model with Modern Architecture
#'
#' @param data joint_data object
#' @param spec JointModelSpec object or formula list
#' @param proteins Character vector of protein names to model
#' @param progress Logical for progress bar
#' @param ... Additional arguments
#'
#' @return joint_model_results object
#' @export
fit_joint_models <- function(
  data,
  spec,
  proteins = NULL,
  progress = TRUE,
  ...
) {
  # Validate inputs
  if (!inherits(data, "joint_data")) {
    cli::cli_abort("data must be a joint_data object")
  }

  if (!inherits(spec, "JointModelSpec")) {
    if (is.list(spec)) {
      # Convert list to spec
      spec <- JointModelSpec$new(
        longitudinal_formula = spec$longitudinal,
        survival_formula = spec$survival,
        association = spec$association %||% "value"
      )
    } else {
      cli::cli_abort("spec must be a JointModelSpec object or formula list")
    }
  }

  # Auto-detect proteins if not specified
  if (is.null(proteins)) {
    # Find numeric columns that aren't covariates
    formula_vars <- all.vars(spec$longitudinal_formula)
    proteins <- setdiff(
      names(dplyr::select(data$longitudinal, where(is.numeric))),
      c(formula_vars, data$id_var, data$time_var)
    )
  }

  cli::cli_alert_info(glue::glue(
    "fit_joint_models: {length(proteins)} protein(s)"
  ))

  if (progress) {
    pb <- progress::progress_bar$new(
      format = "[:bar] :percent :current/:total ETA: :eta",
      total  = length(proteins),
      clear  = FALSE
    )
  }

    results <- furrr::future_map(
      proteins,
      function(protein) {
        if (progress) pb$tick()
        cli::cli_alert_info(glue::glue("Starting protein {protein} at {Sys.time()}"))
        safely_fit_joint_model(data, spec, protein)
      },
      .options = furrr::furrr_options(seed = TRUE)
    )

  cli::cli_alert_success("Finished fit_joint_models")

  joint_model_results(results, proteins, spec, data)
}
