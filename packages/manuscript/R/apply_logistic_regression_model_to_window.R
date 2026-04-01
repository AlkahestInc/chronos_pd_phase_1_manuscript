## SPDX-FileCopyrightText: 2026 Alkahest, Inc. <chronospd_data_access@alkahest.com>
##
## SPDX-License-Identifier: MIT

apply_logistic_regression_model_to_window <- function(.data, model, min_samples = 10, exposure = "ratio") {
  
  # Check if we have enough data
  if (nrow(.data) < min_samples) {
    return(tibble::tibble(
      term = character(),
      estimate = double(),
      std.error = double(),
      statistic = double(),
      p.value = double(),
      conf.low = double(),
      conf.high = double(),
      n_samples = integer(),
      AIC = double(),
      BIC = double(),
      log_or = double(),
      neg_log_p = double(),
      significant = logical()
    ))
  }
  
  model_result <- .data |>
    #####
    # Select only the fields we need for the modelling
    #####
    dplyr::select(
      tidyselect::all_of(all.vars(model))
    ) |>
    #####
    # Drop any data that is incomplete
    #####
    tidyr::drop_na() |>
    #####
    # Apply the logistic regression model to the data
    #####
    apply_logistic_regression_model(.data, model = model, exposure = exposure) |>
    #####
    # Derive some statistics-based columns
    #####
    dplyr::mutate(
      log_or = log2(odds_ratio),
      neg_log_p = -log10(p.value),
      significant = p.value < 0.05
    )
  
  #####
  # Add balance metrics for covariates, as calculated by cobalt, to the results
  #####
  lhs_terms <- all.vars(model[[2]])
  rhs_terms <- all.vars(model[[3]])
  cbm <- cobalt::bal.tab(
    # Only include non-exposure, covariates from the rhs of the model
    reformulate(setdiff(rhs_terms, exposure), lhs_terms),
    data = .data,
    s.d.denom = "pooled",
  )
  
  model_result <- model_result |>
    dplyr::bind_cols(
      cbm$Balance |>
        tibble::as_tibble(rownames = "variable") |>
        dplyr::select(
          variable,
          Diff.Un,
        ) |>
        tidyr::pivot_wider(
          names_from = variable,
          values_from = Diff.Un,
          names_glue = "cobalt_{variable}"
        )
    )
  
  return(model_result)
}
