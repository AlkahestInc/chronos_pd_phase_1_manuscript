## SPDX-FileCopyrightText: 2026 Alkahest, Inc. <chronospd_data_access@alkahest.com>
##
## SPDX-License-Identifier: MIT

#' Refit a mixed model robustly (avoids stats::update call-env issues)
#'
#' @param fit A fitted model (lmer/lmerTest or nlme).
#' @param REML Logical. TRUE for REML, FALSE for ML (lmer) or method="ML" (nlme).
#' @return A refitted model with the requested estimation mode.
refit_mixed <- function(fit, REML = TRUE) {
  if (inherits(fit, "lmerMod")) {
    # Refit using stored formula and the original model frame
    lme4::lmer(formula = stats::formula(fit), data = fit@frame, REML = REML)
  } else if (inherits(fit, "lme")) {
    # nlme supports update() safely
    method <- if (REML) "REML" else "ML"
    stats::update(fit, method = method)
  } else {
    stop("Unsupported fit class in .refit_mixed().")
  }
}
