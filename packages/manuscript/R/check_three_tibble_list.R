## SPDX-FileCopyrightText: 2026 Alkahest, Inc. <chronospd_data_access@alkahest.com>
##
## SPDX-License-Identifier: MIT

#' Check that a list has the required components and basic structures
check_three_tibble_list <- function(three_tibble_list) {
  required_components <- c("prot_data", "meta_data", "prot_annot")

  missing <- setdiff(required_components, names(three_tibble_list))
  if (length(missing) > 0) {
    stop(sprintf("three_tibble_list is missing components: %s",
                 paste(missing, collapse = ", ")))
  }
  for (nm in required_components) {
    obj <- three_tibble_list[[nm]]
    if (!is.data.frame(obj)) {
      stop(sprintf("three_tibble_list$%s must be a data.frame/tibble.", nm))
    }
  }
  invisible(TRUE)
}
