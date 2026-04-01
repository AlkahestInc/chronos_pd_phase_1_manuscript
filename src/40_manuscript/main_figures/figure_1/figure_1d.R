## SPDX-FileCopyrightText: 2026 Alkahest, Inc. <chronospd_data_access@alkahest.com>
##
## SPDX-License-Identifier: MIT

#!/usr/bin/env Rscript

#' Generate Figure 1d: Real-World Data (RWD) Diagnosis Code Timeline
#'
#' This script creates an RWD plot visualization showing the timeline of ICD-10
#' diagnosis codes for Parkinson's Disease patients relative to their estimated
#' disease onset.
#'
#' The plot displays diagnosis codes (NOT samples/donations):
#' - All ICD10 diagnosis codes (background layer, light blue)
#' - G20 (Parkinson's) ICD10 codes (foreground layer, dark red)
#' - Vertical reference line at disease onset (x = 0)
#' - Time window: -15 to +10 years relative to onset
#'
#' Note: This is NOT a tornado plot. Tornado plots show samples; this shows diagnosis codes.
#'
#' Input: Parquet file with diagnosis code data from data publishing figure_1d.R
#' Output: PDF RWD plot showing diagnosis code timeline

here::i_am("src/40_manuscript/main_figures/figure_1/figure_1d.R")

#### Configuration ####
input_dataset_file <- glue::glue(
  "{get_config()$input_dir}/name=chronos_pd-rwd_plot/version=1/dataset.parquet"
)

output_figure_file <- glue::glue(
  "{get_config()$output_dir}/main_figures/figure_1d.pdf"
)
# Create output directory if needed
dir.create(dirname(output_figure_file), showWarnings = FALSE, recursive = TRUE)

#### Load Data ####
logger::log_info("Loading RWD plotting data from {input_dataset_file}", input_dataset_file = input_dataset_file)

# Validate input file exists
if (!file.exists(input_dataset_file)) {
  stop("Input file not found: ", input_dataset_file, "\n",
       "Please run the data generation script first:\n",
       "  Rscript /path/to/data_publishing/src/40_manuscript/main_figures/figure_1/figure_1d.R")
}

rwd_plotting_data <- arrow::read_parquet(input_dataset_file)

logger::log_info("Loaded {nrow(rwd_plotting_data)} diagnosis records for {n_donors} donors",
  nrow = nrow(rwd_plotting_data),
  n_donors = dplyr::n_distinct(rwd_plotting_data$PERSON_ID)
)

#### Prepare Data for Plotting ####
logger::log_info("Filtering data to 2014-2025 date range")

# Filter to relevant time period and select plotting columns
plotting_data <- rwd_plotting_data |>
  # dplyr::filter(
  #   DIAGNOSIS_CODE_DATE >= "2014-01-01",
  #   DIAGNOSIS_CODE_DATE < "2025-01-01"
  # ) |>
  dplyr::select(
    x = DELTA_YEARS,        # Years from disease onset (x-axis)
    y = PERSON_ID,          # Patient identifier (y-axis, already factored with proper order)
    group = INCLUSION_CODE  # G20 vs other ICD10 codes (for layer coloring)
  )

logger::log_info("Filtered to {nrow(plotting_data)} records", nrow = nrow(plotting_data))

# Validate we have data after filtering
if (nrow(plotting_data) == 0) {
  stop("No data remaining after date filtering. Check input data date range.")
}

#### Define Layer Configuration ####
# Two layers:
#   1. Background: All ICD10 codes (light, subtle)
#   2. Foreground: G20 Parkinson's codes (darker, more prominent)
layer_config <- tibble::tribble(
  ~group,              ~size,  ~alpha,  ~colour,      ~shape,
  "ICD10 codes",        0.075,  0.05,   "steelblue",  "circle",
  "G20 ICD10 codes",    0.075,  0.2,    "darkred",    "circle"
)

#### Generate RWD Plot ####
logger::log_info("Generating RWD diagnosis code plot")

rwd_p <- plotting_data |>
  rwd_plot(
    layer_config = layer_config,
    title = "Real-World Data",
    subtitle = "",
    caption = "",
    xlab = "Years Since Estimated Disease Onset",
    ylab = "Donors"
  ) +
  # Constrain x-axis to -15 to +10 years
  ggplot2::coord_cartesian(
    xlim = c(-15L, 10L),
    ylim = c(NA, NA)
  ) +
  # Position legend inside plot area (top-left)
  ggplot2::theme(
    legend.position = "inside",
    legend.position.inside = c(0.15, 0.9)
  )

#### Save Figure ####
logger::log_info("Saving figure to {output_figure_file}", output_figure_file = output_figure_file)

# Save as PDF
ggplot2::ggsave(
  filename = output_figure_file,
  plot = rwd_p,
  width = 8.0,
  height = 4.0,
  units = "in",
  dpi = 600,
  bg = "white"
)

# Save as TIFF
ggplot2::ggsave(
  filename = sub("\\.pdf$", ".tiff", output_figure_file),
  plot = rwd_p,
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
