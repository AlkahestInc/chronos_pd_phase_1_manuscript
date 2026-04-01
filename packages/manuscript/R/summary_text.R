## SPDX-FileCopyrightText: 2026 Alkahest, Inc. <chronospd_data_access@alkahest.com>
##
## SPDX-License-Identifier: MIT

summary_text <- function(data, column, digits = 2) {
  x <- dplyr::pull(data, {{ column }})

  stats <- summary(x)

  stats <- c(
    Min = stats[["Min."]],
    `25th percentile` = stats[["1st Qu."]],
    Median = stats[["Median"]],
    Mean = stats[["Mean"]],
    `75th percentile` = stats[["3rd Qu."]],
    Max = stats[["Max."]]
  )
  stats <- round(stats, digits)

  values <- formatC(stats, format = "f", digits = digits, big.mark = ",")

  name_width  <- max(nchar(names(stats)))
  value_width <- max(nchar(values))

  fmt <- paste0("%-", name_width, "s  %", value_width, "s")

  label_text <- paste(
    sprintf(fmt, names(stats), values),
    collapse = "\n"
  )

  return(label_text)
}
