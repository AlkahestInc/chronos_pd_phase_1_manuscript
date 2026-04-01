## SPDX-FileCopyrightText: 2026 Alkahest, Inc. <chronospd_data_access@alkahest.com>
##
## SPDX-License-Identifier: MIT

#' Split a string on top-level '+' (ignoring '+' inside parentheses)
#' @keywords internal
split_top_level_plus <- function(x) {
  if (length(x) > 1) x <- paste(x, collapse = " + ")
  chars <- strsplit(x, "", fixed = TRUE)[[1]]
  pieces <- character()
  buf <- character()
  depth <- 0L
  for (ch in chars) {
    if (ch == "(") depth <- depth + 1L
    if (ch == ")") depth <- depth - 1L
    if (ch == "+" && depth == 0L) {
      pieces <- c(pieces, paste0(buf, collapse = ""))
      buf <- character()
    } else {
      buf <- c(buf, ch)
    }
  }
  pieces <- c(pieces, paste0(buf, collapse = ""))
  trimws(pieces)
}
