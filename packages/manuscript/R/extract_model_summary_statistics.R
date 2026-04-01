## SPDX-FileCopyrightText: 2026 Alkahest, Inc. <chronospd_data_access@alkahest.com>
##
## SPDX-License-Identifier: MIT

#' extract_model_summary_statistics
#'
#' Extract fixed-effect summary statistics from a fitted model (nlme or lmer/lmerTest).
#'
#' @param fit A fitted model.
#' @return Tibble with columns: coef.names, Estimate, Std.Error, DF, t.value, p.value, cohensD
#' @export
extract_model_summary_statistics <- function(fit) {
  if (inherits(fit, "lme")) {
    ttab <- summary(fit) |>
      coefficients() |>
      tibble::as_tibble(rownames = "coef.names") |>
      dplyr::rename(
        Estimate = "Value",
        Std.Error = "Std.Error",
        DF = "DF",
        t.value = "t-value",
        p.value = "p-value"
      )
  } else if (inherits(fit, "lmerMod")) { # merModLmerTest
    coefs <- summary(fit)$coefficients
    ttab <- tibble::as_tibble(as.data.frame(coefs), rownames = "coef.names") |>
      dplyr::rename(
        Estimate = "Estimate",
        Std.Error = "Std. Error",
        DF = "df",
        t.value = "t value",
        p.value = "Pr(>|t|)"
      )
  } else if (inherits(fit, "merMod")) {
    coefs <- lme4::fixef(fit)
    se <- sqrt(diag(as.matrix(vcov(fit))))
    tvals <- coefs / se
    ttab <- tibble::tibble(
      coef.names = names(coefs),
      Estimate = as.numeric(coefs),
      Std.Error = as.numeric(se),
      DF = NA_real_,
      t.value = as.numeric(tvals),
      p.value = NA_real_
    )
  } else {
    stop("Unsupported fit class in .tidy_fixed_table().")
  }

  dd <- effectsize::t_to_d(ttab$t.value, ttab$DF)
  ttab$cohensD <- dd$d
  ttab
}
