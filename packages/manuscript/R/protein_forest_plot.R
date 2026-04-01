## SPDX-FileCopyrightText: 2026 Alkahest, Inc. <chronospd_data_access@alkahest.com>
##
## SPDX-License-Identifier: MIT

protein_forest_plot <- function(dataset, cohort_order) {
  required_cols <- c("gene_symbol", "lower_ci", "upper_ci", "cohort", "log2FC", "p_val")
  required_cols <- c("cohort", "log2FC")

  if (!all(required_cols %in% names(dataset))) {
    stop(
      "Missing required columns: ",
      paste(setdiff(required_cols, names(dataset)), collapse = ", "),
      call. = FALSE
    )
  }

  stripe_df <- data.frame(
    cohort = factor(cohort_order, levels = cohort_order),
    xmin = seq_along(cohort_order) - 0.5,
    xmax = seq_along(cohort_order) + 0.5,
    fill = rep_len(c("#FFFFFF", "#DDDDDD"), length(cohort_order))
  )

  p <- dataset |>
    dplyr::mutate(
      cohort = factor(cohort, levels = cohort_order)
    ) |>
    ggplot2::ggplot(
      ggplot2::aes(
        x = cohort,
        y = log2FC,
        colour = colour
      )
    ) +
    ggplot2::geom_rect(
      data = stripe_df,
      ggplot2::aes(
        xmin = xmin,
        xmax = xmax,
        ymin = -Inf,
        ymax = Inf,
        fill = fill
      ),
      alpha = 0.5,
      inherit.aes = FALSE
    ) +
    ggplot2::scale_fill_identity() +
    ggplot2::geom_point() +
    ggplot2::geom_errorbar(
      ggplot2::aes(
        ymin = lower_ci,
        ymax = upper_ci
      ),
      width = 0.2
    ) +
    ggplot2::scale_color_manual(
      values = c(
        " " = "#BBBBBB",
        "*" = "#666666",
        "**" = "#333333",
        "***" = "#000000"
      ),
      guide = "none"
    ) +
    ggplot2::scale_x_discrete(drop = FALSE) +
    ggplot2::scale_y_continuous(
      limits = function(x) {
        max_x <- max(x, na.rm = TRUE)
        min_x <- min(x, na.rm = TRUE) * 1.05

        if (max_x > 0) {
          c(min_x, max_x * 1.25) # expand by 15%
        } else {
          c(min_x, max_x)
        }
      }
    ) +
    ggplot2::geom_hline(
      yintercept = 0,
      linetype = "twodash",
      color = "coral2",
      linewidth = 0.5
    ) +
    ggplot2::geom_text(
      ggplot2::aes(
        label = colour,
        y = upper_ci
      ),
      vjust = 0.75,
      hjust = -0.3,
      fontface = "bold"
    ) +
    ggplot2::labs(
      x = "Cohort",
      y = "log2FC"
    ) +
    ggplot2::facet_wrap(~gene_symbol, nrow = 1, scales = "free_x") +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1.0),
      panel.spacing.x = grid::unit(1, "lines")
    ) +
    ggplot2::coord_flip()

  return(p)
}
