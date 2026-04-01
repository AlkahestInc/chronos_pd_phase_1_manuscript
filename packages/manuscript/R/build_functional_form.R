## SPDX-FileCopyrightText: 2026 Alkahest, Inc. <chronospd_data_access@alkahest.com>
##
## SPDX-License-Identifier: MIT

#' Build functional form for association structure
#' @keywords internal
build_functional_form <- function(protein, association) {
  forms <- list()

  if ("value" %in% association) {
    forms$value <- as.formula(paste("~ value(", protein, ")"))
  }
  if ("slope" %in% association) {
    forms$slope <- as.formula(paste("~ slope(", protein, ")"))
  }
  if ("area" %in% association) {
    forms$area <- as.formula(paste("~ area(", protein, ")"))
  }
  if ("shared_random" %in% association) {
    forms$shared <- as.formula(paste("~ shared(", protein, ")"))
  }

  if (length(forms) == 0) {
    # Default to value if nothing specified
    return(as.formula(paste("~ value(", protein, ")")))
  } else if (length(forms) == 1) {
    return(forms[[1]])
  } else {
    # Combine multiple associations
    formula_str <- paste(sapply(forms, deparse), collapse = " + ")
    formula_str <- gsub("~", "", formula_str)
    return(as.formula(paste("~", formula_str)))
  }
}
