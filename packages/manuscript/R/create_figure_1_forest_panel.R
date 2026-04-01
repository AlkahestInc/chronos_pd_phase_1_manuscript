## SPDX-FileCopyrightText: 2026 Alkahest, Inc. <chronospd_data_access@alkahest.com>
##
## SPDX-License-Identifier: MIT

#' Create and save Figure 1 forest plot panel
#'
#' @param chronos_file Path to Chronos-PD Phase 1 results file
#' @param us_file Path to US Established-PD results file
#' @param stage Stage to plot ("before" or "after")
#' @param output_figure_file Path to save the output PDF figure (optional)
#' @param max_codes Maximum number of codes to plot (default: NULL for no limit)
#' @param width Width of the figure in inches (default: 8)
#' @param height Height of the figure in inches (default: 10)
#' @param dpi Resolution in dots per inch (default: 600)
#' @param ... Additional arguments passed to create_forest_plot()
#'
#' @return A list with the plot and plotting data
create_figure_1_forest_panel <- function(
    chronos_file,
    us_file,
    stage,
    output_figure_file = NULL,
    max_codes = NULL,
    width = 8,
    height = 10,
    dpi = 600,
    ...
) {

  # Load and combine data
  data_list <- load_and_combine_pd_data(chronos_file, us_file)

  # Identify significant codes
  significant_codes <- identify_significant_codes_for_stage(
    data_list$res_chronos,
    data_list$res_us,
    stage,
    max_codes
  )

  # Prepare plotting data
  plotting_data <- prepare_forest_plot_data(
    data_list$all_res,
    significant_codes,
    stage
  )

  # Create title based on stage
  title <- dplyr::case_when(
    stage == "before" ~ "Pre-Onset",
    stage == "after" ~ "Post-Onset",
    TRUE ~ stage
  )

  # Create plot
  p <- create_forest_plot(plotting_data, title, ...)

  # Save plot if output file specified
  if (!is.null(output_figure_file)) {
    dir.create(dirname(output_figure_file), showWarnings = FALSE, recursive = TRUE)
    ggplot2::ggsave(
      filename = output_figure_file,
      plot = p,
      width = width,
      height = height,
      units = "in",
      dpi = dpi,
      bg = "white"
    )
  }

  return(p)
}
