## SPDX-FileCopyrightText: 2026 Alkahest, Inc. <chronospd_data_access@alkahest.com>
##
## SPDX-License-Identifier: MIT

plot_volcano <- function(tbl, x, y = "neg_log_p", colour = "significant", title = "Volcano Plot of Logistic Regression Results", x_lab = expression(log[2]("Odds Ratio")), y_lab = expression(-log[10]("p-value"))) {
  required_cols <- c(
    x, y, colour
  )
  missing_cols <- setdiff(required_cols, names(tbl))
  
  if (length(missing_cols) > 0L) {
    stop(
      sprintf("Missing required column(s): %s", paste(missing_cols, collapse = ", ")),
      call. = FALSE
    )
  }
  
  p <- tbl |>
    ggplot2::ggplot(
      ggplot2::aes(
        x = .data[[x]],
        y = .data[[y]],
        colour = .data[[colour]]
      )
    ) +
    ggplot2::geom_point(
      alpha = 0.7,
      shape = 16
    ) +
    ggplot2::geom_vline(
      xintercept = 0,
      linetype = "dashed",
      colour = "grey50"
    ) +
    ggplot2::scale_color_manual(
      values = c("FALSE" = "grey70", "TRUE" = "#D55E00")
    ) +
    ggplot2::labs(
      x = x_lab,
      y = y_lab,
      colour = "Significant",
      title = title
    ) +
    ggplot2::theme_minimal(base_size = 13)
  
  return(p)
}
