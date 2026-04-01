## SPDX-FileCopyrightText: 2026 Alkahest, Inc. <chronospd_data_access@alkahest.com>
##
## SPDX-License-Identifier: MIT

#!/usr/bin/env Rscript

#' Generate Figure 1e: Tornado Plot - Sample Collection Timeline
#'
#' This script creates a tornado plot visualization showing the timeline of plasma sample
#' collections for Parkinson's Disease patients relative to their estimated disease onset.
#'
#' The tornado plot displays sample/donation collection dates (NOT diagnosis codes):
#' - Plasma samples (background layer)
#' - Plasma samples selected for molecular profiling (foreground layer, CHRONOS-PD Phase 1)
#' - Vertical reference line at disease onset (x = 0)
#' - Time window: -15 to +10 years relative to onset
#'
#' This is a TRUE tornado plot showing sample collection events.
#'
#' Input: Parquet file with sample timeline data from data publishing figure_1d.R
#' Output: PDF tornado plot showing sample collection timeline

here::i_am("src/40_manuscript/main_figures/figure_1/figure_1e.R")

#### Configuration ####
input_dataset_file <- glue::glue(
  "{get_config()$input_dir}/name=chronos_pd-tornado_plot/version=1/dataset.parquet"
)

output_figure_file <- glue::glue(
  "{get_config()$output_dir}/main_figures/figure_1e.pdf"
)
dir.create(dirname(output_figure_file), showWarnings = FALSE, recursive = TRUE)

#### Load Data ####
logger::log_info("Loading tornado plot data from {input_dataset_file}",
  input_dataset_file = input_dataset_file
)

if (!file.exists(input_dataset_file)) {
  stop("Input file not found: ", input_dataset_file, "\n",
       "Please run the data generation script first:\n",
       "  Rscript /path/to/data_publishing/src/40_manuscript/main_figures/figure_1/figure_1d.R")
}

tornado_plot_data <- arrow::read_parquet(input_dataset_file)

logger::log_info("Loaded {nrow(tornado_plot_data)} sample records for {n_donors} donors",
  nrow = nrow(tornado_plot_data),
  n_donors = dplyr::n_distinct(tornado_plot_data$PERSON_ID)
)

#### Prepare Data for Plotting ####
logger::log_info("Preparing sample data for tornado plot")

plotting_data <- tornado_plot_data |>
  dplyr::select(
    x = DELTA_YEARS,  # Years from disease onset (x-axis)
    y = PERSON_ID,    # Patient identifier (y-axis, already factored with proper order)
    group = group     # Phase 1 vs regular samples (for layer coloring)
  )

logger::log_info("Prepared {nrow(plotting_data)} sample records", nrow = nrow(plotting_data))

if (nrow(plotting_data) == 0) {
  stop("No sample data available. Check input data.")
}

#### Define Layer Configuration ####
# Two layers:
#   1. Background: Regular plasma samples (light, subtle)
#   2. Foreground: Phase 1 samples selected for profiling (darker, more prominent)
layer_config <- tibble::tribble(
  ~group,                                                   ~size,  ~alpha,  ~colour,      ~shape,
  "Plasma samples",                                         0.25,   0.05,    "#5580B0",    "circle small",
  "Plasma samples selected\nfor molecular profiling",       1.0,    1.0,     "#000085",    "circle small"
)

#### Generate Tornado Plot ####
logger::log_info("Generating tornado plot")

tornado_p <- plotting_data |>
  rwd_plot(
    layer_config = layer_config,
    title = "Biospecimens",
    subtitle = "",
    caption = "",
    xlab = "Years Since Estimated Disease Onset",
    ylab = "Donors"
  ) +
  ggplot2::coord_cartesian(
    xlim = c(-15L, 10L),
    ylim = c(NA, NA)
  ) +
  ggplot2::theme(
    legend.position = "inside",
    legend.position.inside = c(0.81, 0.26)
  )

#### Save Figure ####
logger::log_info("Saving figure to {output_figure_file}", output_figure_file = output_figure_file)

ggplot2::ggsave(
  filename = output_figure_file,
  plot = tornado_p,
  width = 8.0,
  height = 4.0,
  units = "in",
  dpi = 600,
  bg = "white"
)

# Save as TIFF
ggplot2::ggsave(
  filename = sub("\\.pdf$", ".tiff", output_figure_file),
  plot = tornado_p,
  device = "tiff",
  width = 8.0,
  height = 4.0,
  units = "in",
  dpi = 600,
  bg = "white",
  compression = "lzw"
)

logger::log_info("Figure saved successfully")
logger::log_info("Output: {output_figure_file}", output_figure_file = output_figure_file)
logger::log_info("Dimensions: 8.0 x 4.0 inches at 600 DPI")
