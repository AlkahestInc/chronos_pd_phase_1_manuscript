## SPDX-FileCopyrightText: 2026 Alkahest, Inc. <chronospd_data_access@alkahest.com>
##
## SPDX-License-Identifier: MIT

#' explain_by_LME
#'
#' Explanatory Linear Mixed Effects modeling.
#' Supports **random effects as lmer-style strings** (recommended) or legacy structure for nlme.
#'
#' @param protein_id Character. Response column name in `three_tibble_list$prot_data`.
#' @param three_tibble_list See `derive_model_data()`.
#' @param fixed_effects List/character vector of fixed RHS terms (e.g., `"age + sex + x*y"`).
#' @param random_effects Character (lmer-style RHS) **or** legacy list-of-lists (only for `nlme`).
#' @param package One of `"lmerTest"` (default), `"lmer"`, or `"nlme"`.
#' @param ... Passed to the underlying fitter.
#'
#' @return A fitted model (`merModLmerTest`/`merMod`/`nlme::lme`).
#' @export
explain_by_LME <- function(
    protein_id = "",
    three_tibble_list,
    fixed_effects = list(),
    random_effects = NULL,
    package = "lmerTest",
    ...
) {
  model_data <- derive_model_data(protein_id, three_tibble_list, fixed_effects, random_effects)

  if (package %in% c("lmer", "lmerTest")) {
    if (!is.null(random_effects) && !is.character(random_effects)) {
      stop("For package '", package, "', supply `random_effects` as lmer-style character string(s), e.g. '1 | Donor + (0 + slope | Donor) + 1 | Plate'.")
    }
    # Build formula and pass it directly (we won't rely on stats::update later)
    fit <- if (identical(package, "lmerTest")) {
      lmerTest::lmer(formula = build_lmer_formula(protein_id, fixed_effects, random_effects),
                     data = model_data, ...)
    } else {
      lme4::lmer(formula = build_lmer_formula(protein_id, fixed_effects, random_effects),
                 data = model_data, ...)
    }
  } else if (identical(package, "nlme")) {
    if (!is.list(random_effects) || is.character(random_effects)) {
      stop("For package 'nlme', please supply `random_effects` as a legacy list-of-lists as in the original implementation.")
    }
    fixed_effect_formula <- paste0(protein_id, " ~ ",
                                   if (length(fixed_effects) > 0) paste(unlist(fixed_effects), collapse = " + ") else "1")
    random_effect_formula <- paste0("~",
                                    paste(random_effects[["variables"]], collapse = "+"),
                                    "|",
                                    paste(random_effects[["group"]], collapse = "/"))
    fit <- nlme::lme(
      data = model_data,
      fixed = stats::as.formula(fixed_effect_formula),
      random = stats::as.formula(random_effect_formula),
      ...
    )
  } else {
    stop("Invalid package specified. Use 'nlme', 'lmer', or 'lmerTest'.")
  }

  fit
}
