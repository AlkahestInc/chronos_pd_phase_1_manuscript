## SPDX-FileCopyrightText: 2026 Alkahest, Inc. <chronospd_data_access@alkahest.com>
##
## SPDX-License-Identifier: MIT

apply_linear_regression_model_to_data <- function(.data, model, min_samples = 10, exposure = "groupcase.eid") {
  
  # Check if we have enough data
  if (nrow(.data) < min_samples) {
    return(tibble::tibble(
      term = character(),
      beta = double(),
      std.error = double(),
      statistic = double(),
      p.value = double(),
      conf.low = double(),
      conf.high = double(),
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
    apply_linear_regression_model(.data, model = model, exposure = exposure) |>
    #####
    # Derive some statistics-based columns
    #####
    dplyr::mutate(
      neg_log_p = -log10(p.value),
      significant = p.value < 0.05
    )
  
  #####
  # Add balance metrics for covariates, as calculated by cobalt, to the results
  #####
  cbm <- cobalt::bal.tab(
    model,
    data = .data,
    # s.d.denom = "pooled",
  )
  
  model_result <- model_result |>
    dplyr::bind_cols(
      cbm$Balance |>
        tibble::as_tibble(rownames = "variable") |>
        dplyr::select(
          variable,
          Corr.Un,
        ) |>
        tidyr::pivot_wider(
          names_from = variable,
          values_from = Corr.Un,
          names_glue = "cobalt_{variable}"
        )
    )
  
  return(model_result)
}
