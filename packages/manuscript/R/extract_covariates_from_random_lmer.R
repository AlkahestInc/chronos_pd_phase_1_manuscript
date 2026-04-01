## SPDX-FileCopyrightText: 2026 Alkahest, Inc. <chronospd_data_access@alkahest.com>
##
## SPDX-License-Identifier: MIT

#' Extract covariate names (variables & groups) from random-effect lmer-style RHS strings
#' @keywords internal
extract_covariates_from_random_lmer <- function(random_str) {
  if (length(random_str) == 0 || is.null(random_str)) return(character())
  rhs <- paste(random_str, collapse = " + ")
  parts <- split_top_level_plus(rhs)
  vars <- character()
  for (p in parts) {
    p0 <- trimws(p)
    if (startsWith(p0, "(") && endsWith(p0, ")")) p0 <- substr(p0, 2, nchar(p0) - 1)
    sp <- strsplit(p0, "\\|")[[1]]
    if (length(sp) != 2) next
    lhs <- trimws(sp[1])
    grp <- trimws(sp[2])
    lhs_vars <- stringr::str_split(lhs, "\\s*\\+\\s*")[[1]] |>
      trimws() |>
      setdiff(c("0","1"))
    grp_vars <- stringr::str_split(grp, "\\s*/\\s*")[[1]] |>
      trimws()
    vars <- c(vars, lhs_vars, grp_vars)
  }
  unique(vars[nzchar(vars)])
}
