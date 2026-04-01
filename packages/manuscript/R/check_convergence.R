## SPDX-FileCopyrightText: 2026 Alkahest, Inc. <chronospd_data_access@alkahest.com>
##
## SPDX-License-Identifier: MIT

#' Check Convergence for JMbayes2 v0.5-7
#'
#' @description
#' Version-agnostic convergence checking
#'
#' @export
check_convergence <- function(jm_fit, rhat_threshold = 1.05) {
  summary_obj <- summary(jm_fit)

  # Extract Rhat values
  rhat_values <- c()
  if (!is.null(summary_obj$Survival)) {
    rhat_values <- c(rhat_values, summary_obj$Survival[, "Rhat"])
  }
  if (!is.null(summary_obj$Outcome1)) {
    rhat_values <- c(rhat_values, summary_obj$Outcome1[, "Rhat"])
  }

  # Calculate ESS if coda is available and mcmc chains exist
  min_ess <- NA
  ess_info <- list()

  if (!is.null(jm_fit$mcmc) && requireNamespace("coda", quietly = TRUE)) {
    tryCatch(
      {
        # Calculate ESS for different parameter groups
        if (!is.null(jm_fit$mcmc$alphas)) {
          ess_alphas <- coda::effectiveSize(jm_fit$mcmc$alphas)
          ess_info$alphas <- min(ess_alphas)
        }
        if (!is.null(jm_fit$mcmc$betas)) {
          ess_betas <- coda::effectiveSize(jm_fit$mcmc$betas[[1]])
          ess_info$betas <- min(ess_betas)
        }
        if (length(ess_info) > 0) {
          min_ess <- min(unlist(ess_info), na.rm = TRUE)
        }
      },
      error = function(e) {
        # ESS calculation failed
        min_ess <- NA
      }
    )
  }

  list(
    converged = all(rhat_values < rhat_threshold, na.rm = TRUE),
    max_rhat = max(rhat_values, na.rm = TRUE),
    min_ess = min_ess,
    ess_details = ess_info,
    problematic_params = names(rhat_values)[rhat_values >= rhat_threshold]
  )
}
