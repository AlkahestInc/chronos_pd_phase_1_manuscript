## SPDX-FileCopyrightText: 2026 Alkahest, Inc. <chronospd_data_access@alkahest.com>
##
## SPDX-License-Identifier: MIT

#' Create forest plot for comorbidity analysis
#'
#' @param plotting_data Data frame with prepared plotting data
#' @param title Plot title
#' @param x_label Label for x-axis (default: "Log Odds Ratio (95% CI)")
#' @param y_label Label for y-axis (default: "Category")
#' @param color_label Label for color legend (default: "Group")
#' @param point_size Size of points (default: 3)
#' @param point_alpha Point transparency (default: 0.5)
#' @param base_size Base font size (default: 20)
#' @param legend_text_size Legend text size (default: 12)
#' @param color_values Named vector of colors (default: c("blue", "#555555", "MediumSeaGreen"))
#'
#' @return A ggplot2 object
create_forest_plot <- function(
    plotting_data,
    title,
    x_label = "Log Odds Ratio (95% CI)",
    y_label = "Category",
    color_label = "Group",
    point_size = 3,
    point_alpha = 0.5,
    base_size = 20,
    legend_text_size = 12,
    color_values = c("blue", "#555555", "MediumSeaGreen")
) {

  p <- plotting_data |>
    ggplot2::ggplot(
      ggplot2::aes(
        y = forcats::fct_reorder(label_to_plot, dplyr::desc(order_to_plot)),
        x = logOR,
        xmin = LowerCI,
        xmax = UpperCI,
        color = dataset
      )
    ) +
    ggplot2::labs(
      x = x_label,
      y = y_label,
      color = color_label
    ) +
    ggplot2::geom_point(
      size = point_size,
      alpha = point_alpha,
      position = ggplot2::position_dodge(width = 0.5)
    ) +
    ggplot2::scale_color_manual(values = color_values) +
    ggplot2::geom_errorbarh(
      height = 0.15,
      position = ggplot2::position_dodge(width = 0.5)
    ) +
    ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      legend.position = "bottom",
      legend.justification = "left",
      legend.text = ggplot2::element_text(size = legend_text_size)
    ) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
    ggplot2::ggtitle(title)

  return(p)
}
