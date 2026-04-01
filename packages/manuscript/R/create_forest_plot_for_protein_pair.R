## SPDX-FileCopyrightText: 2026 Alkahest, Inc. <chronospd_data_access@alkahest.com>
##
## SPDX-License-Identifier: MIT

# Create forest plot for a specific protein pair
create_forest_plot_for_protein_pair <- function(selected_pair, data) {
  # Filter data for selected pair
  pair_data <- data |>
    dplyr::filter(protein_pair == selected_pair)

  if (nrow(pair_data) == 0) {
    warning(sprintf("No data found for protein pair: %s", selected_pair))
    return(NULL)
  }

  # Get gene symbols for title
  gene1 <- unique(pair_data$gene_symbol1)[1]
  gene2 <- unique(pair_data$gene_symbol2)[1]

  # Convert cohort to factor to preserve order
  # Order: GNPC -> PPMI -> UKB -> Chronos (bottom to top after coord_flip)
  pair_data <- pair_data |>
    dplyr::mutate(
      cohort = factor(cohort, levels = unique(pair_data$cohort)),
      cohort_group = factor(cohort_group, levels = c("GNPC", "PPMI", "UKB", "Chronos"))
    )

  # Calculate positions for background rectangles
  n_cohorts <- length(unique(pair_data$cohort))
  n_gnpc <- sum(pair_data$cohort_group == "GNPC")
  n_ppmi <- sum(pair_data$cohort_group == "PPMI")
  n_ukb <- sum(pair_data$cohort_group == "UKB")
  n_chronos <- sum(pair_data$cohort_group == "Chronos")

  # Calculate y-axis range from the data (with padding for log scale)
  data_range <- range(c(pair_data$or_lower_ci, pair_data$or_upper_ci), na.rm = TRUE)
  y_min <- data_range[1] * 0.5 # Extend range by 50% below
  y_max <- data_range[2] * 1.5 # Extend range by 50% above

  # Create background rectangles data
  # GNPC has 3 independent cohorts (each one line), so create 3 separate boxes with alternating colors
  bg_rects <- data.frame(
    xmin = c(0.5, 1.5, 2.5, n_gnpc + 0.5, n_gnpc + n_ppmi + 0.5, n_gnpc + n_ppmi + n_ukb + 0.5),
    xmax = c(1.5, 2.5, 3.5, n_gnpc + n_ppmi + 0.5, n_gnpc + n_ppmi + n_ukb + 0.5, n_cohorts + 0.5),
    ymin = y_min,
    ymax = y_max,
    group = c("GNPCL", "GNPCJ", "GNPCC", "PPMI", "UKB", "Chronos"),
    fill_color = c("#FFFFFF", "#DDDDDD", "#FFFFFF", "#DDDDDD", "#FFFFFF", "#DDDDDD")
  )

  # Add significance column for coloring
  pair_data <- pair_data |>
    dplyr::mutate(
      significance = pval_to_stars(p_value),
      significance_color = ifelse(significance == "", "ns", significance)
    )

  # Create forest plot
  forest_plot <- pair_data |>
    ggplot2::ggplot(
      ggplot2::aes(
        x = cohort,
        y = odds_ratio,
        color = significance_color
      )
    ) +
    # Background rectangles
    ggplot2::geom_rect(
      data = bg_rects,
      ggplot2::aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = fill_color),
      inherit.aes = FALSE,
      alpha = 0.5,
      show.legend = FALSE
    ) +
    ggplot2::scale_fill_identity() +
    # Points and error bars
    ggplot2::geom_point(
      size = 3,
      position = ggplot2::position_dodge(width = 0.5)
    ) +
    ggplot2::geom_errorbar(
      ggplot2::aes(
        ymin = or_lower_ci,
        ymax = or_upper_ci
      ),
      width = 0.2,
      position = ggplot2::position_dodge(width = 0.5)
    ) +
    # Color scale for significance
    ggplot2::scale_color_manual(
      values = c("ns" = "#BBBBBB", "*" = "#666666", "**" = "#333333", "***" = "#000000"),
      guide = "none"
    ) +
    # Reference line at OR = 1
    ggplot2::geom_hline(
      yintercept = 1,
      linetype = "twodash",
      color = "coral2",
      linewidth = 1
    ) +
    # Significance stars
    ggplot2::geom_text(
      ggplot2::aes(
        label = pval_to_stars(p_value),
        y = or_upper_ci
      ),
      vjust = 0.7,
      hjust = -0.5,
      size = 4,
      fontface = "bold",
      position = ggplot2::position_dodge(width = 0.5),
      show.legend = FALSE
    ) +
    # Flip coordinates and use log scale
    ggplot2::coord_flip() +
    ggplot2::scale_y_log10() +
    # Theme and labels
    ggplot2::theme_minimal() +
    ggplot2::labs(
      title = paste(gene1, gene2, sep = " / "),
      x = NULL,
      y = "Odds Ratio"
    ) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 12),
      axis.text = ggplot2::element_text(size = 9),
      axis.title = ggplot2::element_text(size = 10)
    )

  return(list(
    plot = forest_plot,
    gene1 = gene1,
    gene2 = gene2
  ))
}
