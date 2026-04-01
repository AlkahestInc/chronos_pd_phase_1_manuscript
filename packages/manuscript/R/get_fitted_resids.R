## SPDX-FileCopyrightText: 2026 Alkahest, Inc. <chronospd_data_access@alkahest.com>
##
## SPDX-License-Identifier: MIT

#' Compute fitted values and residuals for lmer or nlme fits
#' Returns a list: response, yhat.mar, yhat.cond, resid.mar, resid.cond
get_fitted_resids <- function(fit, response_name = NULL) {
  if (inherits(fit, "lme")) {
    response <- nlme::getResponse(fit)
    yhat.mar <- stats::fitted(fit, level = 0)
    yhat.cond <- stats::fitted(fit, level = 1)
    resid.mar <- stats::residuals(fit, type = "response", level = 0)
    resid.cond <- stats::residuals(fit, type = "response", level = 1)
  } else if (inherits(fit, "merMod")) {
    frm <- stats::formula(fit)
    resp <- if (is.null(response_name)) as.character(frm[[2]]) else response_name
    yhat.cond <- stats::fitted(fit)               # includes random effects
    yhat.mar <- stats::predict(fit, re.form = NA) # fixed-effects only
    response <- fit@frame[[resp]]
    resid.cond <- stats::residuals(fit, type = "response")
    resid.mar <- response - yhat.mar
  } else {
    stop("Unsupported fit class in .get_fitted_resids().")
  }
  list(
    response = as.numeric(response),
    yhat.mar = as.numeric(yhat.mar),
    yhat.cond = as.numeric(yhat.cond),
    resid.mar = as.numeric(resid.mar),
    resid.cond = as.numeric(resid.cond)
  )
}
