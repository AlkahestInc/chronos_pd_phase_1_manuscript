## SPDX-FileCopyrightText: 2026 Alkahest, Inc. <chronospd_data_access@alkahest.com>
##
## SPDX-License-Identifier: MIT

apply_linear_regression_model <- function(.data, model, exposure, .id = NULL) {
  lhs_terms <- all.vars(stats::formula(model)[[2]])
  rhs_terms <- all.vars(stats::formula(model)[[3]])
  
  required_cols <- c(
    lhs_terms,
    rhs_terms
  )
  missing_cols <- setdiff(required_cols, names(.data))
  
  if (length(missing_cols) > 0L) {
    stop(
      sprintf("Missing required column(s): %s", paste(missing_cols, collapse = ", ")),
      call. = FALSE
    )
  }
  
  n_samples <- nrow(.data)
  
  if (n_samples == 0) {
    result <- tibble::tibble(
      term = character(),
      estimate = double(),
      std.error = double(),
      statistic = double(),
      p.value = double(),
      conf.low = double(),
      conf.high = double(),
      n_samples = integer(),
    )
    
    return(result)
  }
  
  # Fit adjusted linear regression model using tidymodels
  linear_model <- parsnip::linear_reg(mode = "regression") |>
    parsnip::set_engine("glm", family = gaussian()) |>
    parsnip::fit(model, data = .data)
  
  # Extract coefficient for the exposure variable, compute odds ratios and confidence intervals
  linear_model <- broom::tidy(linear_model, conf.int = TRUE) |>
    dplyr::filter(
      term == exposure
    ) |>
    dplyr::rename(
      beta = estimate
    )
  
  # Build return object
  result <- linear_model |>
    dplyr::mutate(n_samples = n_samples)
  
  return(result)
}
