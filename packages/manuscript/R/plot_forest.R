## SPDX-FileCopyrightText: 2026 Alkahest, Inc. <chronospd_data_access@alkahest.com>
##
## SPDX-License-Identifier: MIT

#' Create Modern Forest Plot
#'
#' @param results joint_model_results object
#' @param term Which model term to plot
#' @param n_top Number of top results to show
#' @param show_risk Show risk percentage instead of HR
#' @export
plot_forest <- function(results, term = "value", n_top = 20, show_risk = FALSE) {
  plot_data <- results$summaries |>
    dplyr::filter(term == !!term) |>
    dplyr::arrange(p_value) |>
    dplyr::slice_head(n = n_top)

  if (show_risk) {
    plot_data <- plot_data |>
      dplyr::mutate(
        effect = dplyr::if_else(estimate > 0,
                                (exp(estimate * 0.1) - 1) * 100,
                                (1 - exp(estimate * 0.1)) * 100
        ),
        lower = dplyr::if_else(estimate > 0,
                               (exp(lower_95 * 0.1) - 1) * 100,
                               (1 - exp(lower_95 * 0.1)) * 100
        ),
        upper = dplyr::if_else(estimate > 0,
                               (exp(upper_95 * 0.1) - 1) * 100,
                               (1 - exp(upper_95 * 0.1)) * 100
        ),
        protein = forcats::fct_reorder(protein, effect)
      )

    null_value <- 0
    x_label <- "Risk change (%) per 0.1 unit increase"
  } else {
    plot_data <- plot_data |>
      dplyr::mutate(
        effect = exp(estimate),
        lower = exp(lower_95),
        upper = exp(upper_95),
        protein = forcats::fct_reorder(protein, effect)
      )

    null_value <- 1
    x_label <- "Hazard Ratio"
  }

  ggplot2::ggplot(plot_data, ggplot2::aes(y = protein)) +
    ggplot2::geom_point(ggplot2::aes(x = effect), size = 3) +
    ggplot2::geom_errorbarh(ggplot2::aes(xmin = lower, xmax = upper),
                            height = 0.2
    ) +
    ggplot2::geom_vline(
      xintercept = null_value,
      linetype = "dashed", alpha = 0.5
    ) +
    ggplot2::scale_x_continuous(
      trans = if (show_risk) "identity" else "log10"
    ) +
    ggplot2::labs(
      x = x_label, y = NULL,
      title = glue::glue("Joint Model Results: {term} associations")
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(panel.grid.major.y = ggplot2::element_blank())
}
