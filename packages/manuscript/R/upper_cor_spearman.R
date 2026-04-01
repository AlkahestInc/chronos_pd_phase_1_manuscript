## SPDX-FileCopyrightText: 2026 Alkahest, Inc. <chronospd_data_access@alkahest.com>
##
## SPDX-License-Identifier: MIT

#' Custom upper panel with Spearman correlation for ggpairs plots
#'
#' @param data A data frame
#' @param mapping Aesthetic mapping
#' @param ... Additional arguments
#'
#' @return A ggplot object with correlation text
#' @export
upper_cor_spearman <- function(data, mapping, ...) {

  # Evaluate the mapped columns (mapping$x and mapping$y are quosures)
  x <- rlang::eval_tidy(mapping$x, data = data)
  y <- rlang::eval_tidy(mapping$y, data = data)

  corr <- stats::cor(x, y, method = "spearman", use = "complete.obs")

  # Round first to 3 decimal places
  corr_rounded <- round(corr, 3)

  stars <- symnum(abs(corr_rounded),
                  cutpoints = c(0, 0.3, 0.5, 0.7, 1),
                  symbols = c("", "*", "**", "***"))

  # Dynamic color based on rounded value: red if strong enough, else black
  txt_color <- if (abs(corr_rounded) >= 0.3) "red" else "black"

  GGally::ggally_text(
    label = paste0("rho = ", corr_rounded, "\n", stars),
    mapping = ggplot2::aes(),
    color = txt_color,
    size = 4
  ) +
    ggplot2::theme_void()

}
