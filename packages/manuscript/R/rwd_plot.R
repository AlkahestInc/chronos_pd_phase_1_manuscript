## SPDX-FileCopyrightText: 2026 Alkahest, Inc. <chronospd_data_access@alkahest.com>
##
## SPDX-License-Identifier: MIT

#' Create a Real World Data (RWD) or Tornado Plot
#'
#' Creates a scatter plot for longitudinal event data with customizable layers, colors,
#' and styling. This flexible function can visualize two types of timelines relative to
#' disease onset:
#'
#' 1. **RWD Plot**: ICD-10 diagnosis code occurrences (e.g., Figure 1d)
#' 2. **Tornado Plot**: Sample/donation collection dates (e.g., Figure 1e)
#'
#' The function provides flexible configuration options for different data groups,
#' allowing you to highlight specific events with different colors, sizes, and transparency.
#'
#' @param tbl A data frame containing the data to plot. Must have columns:
#'   \itemize{
#'     \item \code{x}: Numeric time points (e.g., years since disease onset)
#'     \item \code{y}: Factor or character for participant identifiers (y-axis ordering)
#'   }
#'   If \code{layer_config} is provided, must also have a \code{group} column for grouping.
#' @param title Character string. The main title for the plot. Default is "Real-World Data".
#' @param subtitle Character string. The subtitle for the plot. Default is empty string.
#' @param caption Character string. The caption for the plot. Default is empty string.
#' @param xlab Character string. The x-axis label. Default is "Years Since Estimated Disease Onset".
#' @param ylab Character string. The y-axis label. Default is "Donors". The final label will
#'   include the number of unique participants.
#' @param breaks Numeric. The interval for y-axis tick breaks and labels. Default is 500.
#'   Used to avoid overcrowding on plots with many participants.
#' @param pointsize Numeric. The size of points when using default layer configuration. Default is 0.25.
#' @param alpha Numeric. The transparency of points when using default layer configuration.
#'   Value between 0 and 1. Default is 0.025.
#' @param layer_config A data frame specifying custom layer configurations. Must contain columns:
#'   \itemize{
#'     \item \code{group}: Character, matching values in \code{tbl$group}
#'     \item \code{size}: Numeric, point size
#'     \item \code{alpha}: Numeric, transparency (0-1)
#'     \item \code{colour}: Character, color name or hex code
#'     \item \code{shape}: Character, ggplot2 shape name (e.g., "circle", "triangle")
#'   }
#'   If NULL, uses default configuration with single layer.
#'   When provided, \code{tbl} must have a \code{group} column matching the groups in this configuration.
#'
#' @return A ggplot2 object representing the real-world data plot. Can be further customized
#'   with additional ggplot2 layers or saved with \code{ggsave()}.
#'
#' @details
#' The function creates a scatter plot with the following features:
#' \itemize{
#'   \item Customizable point layers based on data groups (useful for highlighting specific diagnosis codes)
#'   \item Automatic y-axis labeling with sample size (n = X,XXX format)
#'   \item Vertical reference line at x = 0 (typically disease onset)
#'   \item Custom theme with white background and minimal styling
#'   \item Legend with enhanced visibility for small/transparent points
#'   \item Discrete y-axis with custom breaks to handle large participant counts
#' }
#'
#' ## Plot Design & Use Cases
#' This function creates timeline visualizations with the following structure:
#' \itemize{
#'   \item Time relative to disease onset (x = 0)
#'   \item One row per participant (y-axis)
#'   \item Points representing events
#'   \item Multiple layers to distinguish event types
#' }
#'
#' **Use Case 1 - RWD Plot (Diagnosis Codes)**:
#' \itemize{
#'   \item Points = ICD-10 diagnosis code occurrences
#'   \item Layers = Different diagnosis types (e.g., G20 vs other codes)
#'   \item Example: Figure 1d
#' }
#'
#' **Use Case 2 - Tornado Plot (Samples)**:
#' \itemize{
#'   \item Points = Sample/donation collection dates
#'   \item Layers = Sample types (e.g., Phase 1 profiled vs regular samples)
#'   \item Example: Figure 1e
#' }
#'
#' The y-axis is treated as discrete with custom breaks and labels. The x-axis shows
#' continuous values with breaks every 5 units from -100 to 100.
#'
#' ## Layer Configuration
#' When using multiple layers (via \code{layer_config}), points are drawn in order,
#' so later groups appear on top. Order groups strategically:
#' \itemize{
#'   \item Background/common diagnoses first (lower alpha, "steelblue")
#'   \item Important/inclusion diagnoses last (higher alpha, "darkred")
#' }
#'
#' @examples
#' \dontrun{
#' # Basic usage with default settings
#' sample_data <- data.frame(
#'   x = rnorm(1000, 0, 10),
#'   y = factor(sample(1:100, 1000, replace = TRUE))
#' )
#'
#' # Basic plot
#' rwd_plot(sample_data)
#'
#' # With custom labels
#' rwd_plot(
#'   sample_data,
#'   title = "Disease Progression Data",
#'   subtitle = "Longitudinal Study Results",
#'   caption = "Source: Clinical Database"
#' )
#'
#' # Example 1: RWD Plot (diagnosis codes)
#' diagnosis_data <- sample_data |>
#'   dplyr::mutate(
#'     group = sample(c("ICD10 codes", "G20 ICD10 codes"), 1000, replace = TRUE)
#'   )
#'
#' rwd_layers <- tibble::tribble(
#'   ~group,              ~size,  ~alpha,  ~colour,      ~shape,
#'   "ICD10 codes",        0.075,  0.05,   "steelblue",  "circle",
#'   "G20 ICD10 codes",    0.075,  0.2,    "darkred",    "circle",
#' )
#'
#' p1 <- rwd_plot(
#'   diagnosis_data,
#'   layer_config = rwd_layers,
#'   title = "Real-World Diagnosis Code Timeline"
#' )
#'
#' # Example 2: Tornado Plot (sample collections)
#' sample_collection_data <- sample_data |>
#'   dplyr::mutate(
#'     group = sample(
#'       c("Plasma samples", "Plasma samples selected\nfor molecular profiling"),
#'       1000,
#'       replace = TRUE
#'     )
#'   )
#'
#' tornado_layers <- tibble::tribble(
#'   ~group,                                                  ~size,  ~alpha,  ~colour,      ~shape,
#'   "Plasma samples",                                         0.5,    0.3,    "steelblue",  "circle",
#'   "Plasma samples selected\nfor molecular profiling",      0.5,    0.8,    "darkred",    "circle",
#' )
#'
#' p2 <- rwd_plot(
#'   sample_collection_data,
#'   layer_config = tornado_layers,
#'   title = "Sample Collection Timeline"
#' )
#'
#' # Further customize the plot
#' p <- p +
#'   ggplot2::coord_cartesian(xlim = c(-15, 10)) +
#'   ggplot2::theme(legend.position = "inside", legend.position.inside = c(0.15, 0.9))
#'
#' # Save to file
#' ggplot2::ggsave("rwd_plot.pdf", plot = p, width = 8, height = 4, dpi = 600)
#' }
#'
#' @importFrom dplyr mutate distinct tally pull
#' @importFrom ggplot2 ggplot aes geom_point scale_color_manual labs scale_y_discrete
#'   theme element_rect element_line element_blank element_text scale_x_continuous
#'   guides guide_legend geom_vline
#' @importFrom glue glue
#' @importFrom tibble tibble
#' @importFrom logger log_fatal
#'
#' @export
rwd_plot <- function(
    tbl,
    title = "Real-World Data",
    subtitle = "",
    caption = "",
    xlab = "Years Since Estimated Disease Onset",
    ylab = "Donors",
    breaks = 500,
    pointsize = 0.25,
    alpha = 0.025,
    layer_config = NULL) {
  # Input validation
  stopifnot(
    "tbl must be a data.frame" = is.data.frame(tbl),
    "title must be a character string" = is.character(title) && length(title) == 1,
    "subtitle must be a character string" = is.character(subtitle) && length(subtitle) == 1,
    "caption must be a character string" = is.character(caption) && length(caption) == 1,
    "xlab must be a character string" = is.character(xlab) && length(xlab) == 1,
    "ylab must be a character string" = is.character(ylab) && length(ylab) == 1,
    "breaks must be a positive number" = is.numeric(breaks) && breaks > 0,
    "pointsize must be a positive number" = is.numeric(pointsize) && pointsize > 0,
    "alpha must be between 0 and 1" = is.numeric(alpha) && alpha >= 0 && alpha <= 1
  )

  # Validate required columns
  required_cols <- c("x", "y")
  missing_cols <- setdiff(required_cols, names(tbl))
  if (length(missing_cols) > 0) {
    stop("tbl missing required columns: ", paste(missing_cols, collapse = ", "))
  }

  # Validate x is numeric
  if (!is.numeric(tbl$x)) {
    stop("Column 'x' must be numeric")
  }

  # Default layer configuration if none provided
  if (is.null(layer_config)) {
    layer_config <- tibble::tibble(
      group = "Other Code",
      size = pointsize,
      alpha = alpha,
      colour = "steelblue",
      shape = "circle"
    )
  } else {
    # Validate layer_config structure
    required_config_cols <- c("group", "size", "alpha", "colour", "shape")
    missing_config_cols <- setdiff(required_config_cols, names(layer_config))
    if (length(missing_config_cols) > 0) {
      stop("layer_config missing required columns: ", paste(missing_config_cols, collapse = ", "))
    }

    # layer config provided, we need a "group" column in tbl
    if (!"group" %in% colnames(tbl)) {
      logger::log_fatal(
        "When providing a 'layer_config', you also need to provide a 'group' column for tbl."
      )
    }

    # Validate that all groups in tbl exist in layer_config
    tbl_groups <- unique(tbl$group)
    config_groups <- layer_config$group
    missing_groups <- setdiff(tbl_groups, config_groups)
    if (length(missing_groups) > 0) {
      stop(
        "tbl contains groups not defined in layer_config: ",
        paste(missing_groups, collapse = ", ")
      )
    }
  }

  # Add default group column if not present
  if (!"group" %in% colnames(tbl)) {
    tbl <- tbl |>
      dplyr::mutate(group = "Other Code")
  }

  # Set the y-axis label to include the number of individuals
  ylab <- glue::glue(
    "{ylab}\nn = {y_count}",
    y_count = tbl |>
      dplyr::distinct(y) |>
      dplyr::tally() |>
      dplyr::pull() |>
      format(big.mark = ",")
  )

  # Initialize plot with invisible base layer (for proper legend)
  p <- tbl |>
    ggplot2::ggplot(
      ggplot2::aes(
        x = x,
        y = y,
        colour = group
      )
    ) +
    ggplot2::geom_point(size = 0.0, alpha = 0.0)

  # Add layers dynamically for each group
  for (i in seq_len(nrow(layer_config))) {
    config <- layer_config[i, ]
    p <- p +
      ggplot2::geom_point(
        data = subset(tbl, group == config$group),
        size = config$size,
        alpha = config$alpha,
        colour = config$colour,
        shape = config$shape
      )
  }

  # Add color scale and labels
  p <- p +
    ggplot2::scale_color_manual(
      values = setNames(layer_config[["colour"]], layer_config[["group"]])
    ) +
    ggplot2::labs(
      title = title,
      subtitle = subtitle,
      caption = caption,
      x = xlab,
      y = ylab
    ) +
    ggplot2::scale_y_discrete(
      breaks = levels(tbl$y)[seq(1L, length(levels(tbl$y)), by = breaks)],
      labels = seq(1L, length(levels(tbl$y)), by = breaks) - 1L
    )

  # Apply custom theme
  p <- p +
    ggplot2::theme(
      panel.background = ggplot2::element_rect(fill = "white"),
      axis.line = ggplot2::element_line(colour = "black")
    ) +
    ggplot2::theme(
      aspect.ratio = 0.5,
      axis.text.y = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      legend.position = "right",
      legend.title = ggplot2::element_blank(),
      plot.caption = ggplot2::element_text(hjust = 0.0),
      plot.caption.position = "plot",
      plot.title = ggplot2::element_text(hjust = 0.5),
      plot.subtitle = ggplot2::element_text(hjust = 0.5)
    ) +
    # Set the x-axis tick interval to every 5 years
    ggplot2::scale_x_continuous(
      breaks = seq(-100L, 100L, by = 5L)
    ) +
    # Enhance legend visibility for small/transparent points
    ggplot2::guides(
      colour = ggplot2::guide_legend(
        override.aes = list(
          alpha = 1.0,
          size = 3.0
        )
      )
    ) +
    # Add vertical reference line at x = 0 (disease onset)
    ggplot2::geom_vline(
      xintercept = 0.0,
      linetype = "solid",
      color = "darkred",
      linewidth = 0.3
    )

  p
}
