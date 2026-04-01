## SPDX-FileCopyrightText: 2026 Alkahest, Inc. <chronospd_data_access@alkahest.com>
##
## SPDX-License-Identifier: MIT

#' Extract Model Summary for JMbayes2 v0.5-7
#'
#' @description
#' Works with the actual structure in your version
#'
#' @export
extract_model_summary <- function(jm_fit, protein) {
  summary_obj <- summary(jm_fit)

  # Extract survival coefficients
  surv_coef <- summary_obj$Survival

  # Find association terms (value, slope, area)
  association_patterns <- paste0(c("value", "slope", "area"), "\\(", protein, "\\)")
  association_rows <- grep(
    paste(association_patterns, collapse = "|"),
    rownames(surv_coef)
  )

  if (length(association_rows) > 0) {
    results <- tibble::tibble(
      protein = protein,
      term = gsub(
        paste0("\\(", protein, "\\)"), "",
        rownames(surv_coef)[association_rows]
      ),
      estimate = surv_coef[association_rows, "Mean"],
      std_error = surv_coef[association_rows, "StDev"],
      lower_95 = surv_coef[association_rows, "2.5%"],
      upper_95 = surv_coef[association_rows, "97.5%"],
      p_value = surv_coef[association_rows, "P"],
      rhat = surv_coef[association_rows, "Rhat"]
    )
  } else {
    # Return NA row if no associations found
    results <- tibble::tibble(
      protein = protein,
      term = "none_found",
      estimate = NA_real_,
      std_error = NA_real_,
      lower_95 = NA_real_,
      upper_95 = NA_real_,
      p_value = NA_real_,
      rhat = NA_real_
    )
  }

  # Add fit statistics
  if (!is.null(jm_fit$fit_stats)) {
    if (!is.null(jm_fit$fit_stats$conditional)) {
      results$dic <- jm_fit$fit_stats$conditional$DIC %||% NA_real_
      results$waic <- jm_fit$fit_stats$conditional$WAIC %||% NA_real_
      results$lpml <- jm_fit$fit_stats$conditional$LPML %||% NA_real_
    } else if (!is.null(jm_fit$fit_stats$marginal)) {
      results$dic <- jm_fit$fit_stats$marginal$DIC %||% NA_real_
      results$waic <- jm_fit$fit_stats$marginal$WAIC %||% NA_real_
      results$lpml <- jm_fit$fit_stats$marginal$LPML %||% NA_real_
    } else {
      results$dic <- NA_real_
      results$waic <- NA_real_
      results$lpml <- NA_real_
    }
  } else {
    results$dic <- NA_real_
    results$waic <- NA_real_
    results$lpml <- NA_real_
  }

  return(results)
}
