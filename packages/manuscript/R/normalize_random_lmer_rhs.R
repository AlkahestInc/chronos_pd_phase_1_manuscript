## SPDX-FileCopyrightText: 2026 Alkahest, Inc. <chronospd_data_access@alkahest.com>
##
## SPDX-License-Identifier: MIT

#' Normalize random-effect character RHS to lmer terms like '(1 | grp)'
#' Accepts inputs like '1 | Donor + 1 | Plate' or '(0 + slope | Donor) + 1 | Plate'
#' @keywords internal
normalize_random_lmer_rhs <- function(random_str) {
  if (length(random_str) > 1) random_str <- paste(random_str, collapse = " + ")
  terms <- split_top_level_plus(random_str)
  terms <- vapply(
    terms,
    function(ti) {
      ti_trim <- trimws(ti)
      if (startsWith(ti_trim, "(") && endsWith(ti_trim, ")")) {
        ti_trim
      } else {
        paste0("(", ti_trim, ")")
      }
    },
    character(1)
  )
  paste(terms, collapse = " + ")
}
