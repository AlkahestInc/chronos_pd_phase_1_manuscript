## SPDX-FileCopyrightText: 2026 Alkahest, Inc. <chronospd_data_access@alkahest.com>
##
## SPDX-License-Identifier: MIT

# ------------------------------------------------------------------------------
# Function: multi_rank
# Description: Rank data frame rows using multiple tie-breaker columns
# Parameters:
#   ... - vectors to rank by (in order of priority)
#   ties.method - method for handling ties (passed to data.table::frank)
# Returns: integer vector of ranks
# ------------------------------------------------------------------------------
multi_rank <- function(..., ties.method = "min") {
  args <- list(...)
  stopifnot(length(unique(sapply(args, length))) == 1)
  
  # Convert to data.table
  dt <- data.table::as.data.table(args)
  
  # Use frank with multiple ranking keys
  return(data.table::frank(dt, ties.method = ties.method))
  
  if (DEV_MODE) {
    # Usage example
    example_df <- data.frame(
      p = c(0, 0, 0, 0.1),
      est = c(10, 3, 10, 100),
      bubu = c(0, 10, 0.1, 100)
    )
    example_df$rank <- multi_rank(example_df$p, -example_df$est, -example_df$bubu)
    example_df <- example_df[order(example_df$rank), ]
  }
}
