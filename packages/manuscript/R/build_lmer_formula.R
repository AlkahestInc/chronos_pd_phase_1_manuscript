## SPDX-FileCopyrightText: 2026 Alkahest, Inc. <chronospd_data_access@alkahest.com>
##
## SPDX-License-Identifier: MIT

#' Build full model formula for lmer/lmerTest
#' @keywords internal
build_lmer_formula <- function(response, fixed_effects, random_str) {
  fixed_rhs <- if (length(fixed_effects) > 0) {
    if (length(fixed_effects) == 1) {
      as.character(fixed_effects)
    } else {
      paste(unlist(fixed_effects, use.names = FALSE), collapse = " + ")
    }
  } else {
    "1"
  }
  rand_rhs <- if (is.null(random_str) || length(random_str) == 0) {
    ""
  } else {
    paste0(" + ", normalize_random_lmer_rhs(random_str))
  }
  stats::as.formula(paste0(response, " ~ ", fixed_rhs, rand_rhs))
}
