## SPDX-FileCopyrightText: 2026 Alkahest, Inc. <chronospd_data_access@alkahest.com>
##
## SPDX-License-Identifier: MIT

#' Custom lower panel with smoothing line for ggpairs plots
#'
#' @param data A data frame
#' @param mapping Aesthetic mapping
#' @param ... Additional arguments
#'
#' @return A ggplot object
#' @export
lower_smooth <- function(data, mapping, ...) {
  ggplot2::ggplot(data = data, mapping = mapping) +
    ggplot2::geom_point(alpha = 0.5, size = 0.5) +
    ggplot2::geom_smooth(method = "loess", se = FALSE, color = "red", size = 0.5, formula = y ~ x) +
    ggplot2::theme_minimal()
}
