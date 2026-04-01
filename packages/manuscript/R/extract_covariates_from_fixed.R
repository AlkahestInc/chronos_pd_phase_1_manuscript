## SPDX-FileCopyrightText: 2026 Alkahest, Inc. <chronospd_data_access@alkahest.com>
##
## SPDX-License-Identifier: MIT

#' Extract covariate names from a vector/list/character of fixed-effect RHS strings
#' We exclude some known function names commonly used in formulas.
#' @keywords internal
extract_covariates_from_fixed <- function(fixed_effects) {
  if (length(fixed_effects) == 0) return(character())
  fxs <- unlist(fixed_effects, use.names = FALSE)
  toks <- stringr::str_split(fxs, pattern = "[^A-Za-z0-9_.]+", simplify = FALSE) |>
    unlist(use.names = FALSE) |>
    unique()
  toks <- toks[nzchar(toks)]
  exclude_fns <- c("I","poly","ns","bs","s","offset","scale",
                   "log","log1p","sqrt","as.factor","factor","cut")
  setdiff(toks, c("0","1", exclude_fns))
}
