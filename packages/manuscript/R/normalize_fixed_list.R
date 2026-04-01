## SPDX-FileCopyrightText: 2026 Alkahest, Inc. <chronospd_data_access@alkahest.com>
##
## SPDX-License-Identifier: MIT

#' Normalize fixed-effect specifications to a character vector of RHS strings
#'
#' Accepts either:
#' - a character vector where each element is a RHS (e.g., "age + sex"),
#' - a list where each element is a character vector of tokens to be joined with " + ".

normalize_fixed_list <- function(fixed_effects_list) {
  if (is.character(fixed_effects_list)) {
    fixed_effects_list
  } else {
    purrr::map_chr(
      fixed_effects_list,
      ~ paste(unlist(.x, use.names = FALSE), collapse = " + ")
    )
  }
}
