## SPDX-FileCopyrightText: 2026 Alkahest, Inc. <chronospd_data_access@alkahest.com>
##
## SPDX-License-Identifier: MIT

apply_logistic_regression_model <- function(.data, model, exposure, .id = NULL) {
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
      AIC = double(),
      BIC = double(),
    )
    
    return(result)
  }
  
  # Fit logistic regression model using tidymodels
  # TODO: look at using broom::glance() instead: https://rpubs.com/DragonflyStats/broom_R_package
  logistic_model <- parsnip::logistic_reg(mode = "classification") |>
    parsnip::set_engine("glm") |>
    parsnip::fit(model, data = .data)
  
  # Extract coefficient for the exposure variable, compute odds ratios and confidence intervals
  odds_ratios <- broom::tidy(logistic_model, exponentiate = TRUE, conf.int = TRUE) |>
    dplyr::filter(
      term == exposure
    ) |>
    dplyr::rename(
      odds_ratio = estimate
    )
  
  # Build return object
  result <- odds_ratios |>
    dplyr::mutate(
      n_samples = n_samples,
      AIC = AIC(logistic_model$fit),
      BIC = BIC(logistic_model$fit),
    )
  
  return(result)
}
