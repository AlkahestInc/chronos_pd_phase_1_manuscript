## SPDX-FileCopyrightText: 2026 Alkahest, Inc. <chronospd_data_access@alkahest.com>
##
## SPDX-License-Identifier: MIT

#' extract_model_diagnostic_statistics
#'
#' Calculate diagnostic statistics for conditional and marginal residuals.
#'
#' @param fit A fitted model (nlme or lmer/lmerTest).
#' @return Tibble with summary diagnostics.
#' @export
extract_model_diagnostic_statistics <- function(fit) {
  resp_name <- if (inherits(fit, "lme")) {
    as.character(fit$terms[[2]])
  } else {
    as.character(stats::formula(fit)[[2]])
  }
  fr <- get_fitted_resids(fit, response_name = resp_name)

  IQRc <- stats::quantile(fr$resid.cond, 0.75) - stats::quantile(fr$resid.cond, 0.25)
  IQRm <- stats::quantile(fr$resid.mar, 0.75) - stats::quantile(fr$resid.mar, 0.25)

  tibble::tibble(
    protein_ID = resp_name,
    mean_protein = mean(fr$response),
    mean_yhat = mean(fr$yhat.cond),
    median_protein = stats::median(fr$response),
    median_yhat = stats::median(fr$yhat.cond),

    shapirotest.pvalue.cond = stats::shapiro.test(fr$resid.cond)$p.value,
    resid.cond.mean = mean(fr$resid.cond),
    resid.cond.corr = stats::cor(fr$yhat.cond, fr$resid.cond, method = "spearman"),
    yhat.cond.corr = stats::cor(fr$response, fr$yhat.cond, method = "spearman"),
    outliers.cond = sum(fr$resid.cond < (mean(fr$resid.cond) - 1.5 * IQRc) |
                          fr$resid.cond > (mean(fr$resid.cond) + 1.5 * IQRc)),

    shapirotest.pvalue = stats::shapiro.test(fr$resid.mar)$p.value,
    resid.mar.mean = mean(fr$resid.mar),
    resid.mar.corr = stats::cor(fr$yhat.mar, fr$resid.mar, method = "spearman"),
    yhat.mar.corr = stats::cor(fr$response, fr$yhat.mar, method = "spearman"),
    outliers = sum(fr$resid.mar < (mean(fr$resid.mar) - 1.5 * IQRm) |
                     fr$resid.mar > (mean(fr$resid.mar) + 1.5 * IQRm))
  )
}
