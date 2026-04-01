## SPDX-FileCopyrightText: 2026 Alkahest, Inc. <chronospd_data_access@alkahest.com>
##
## SPDX-License-Identifier: MIT

here::i_am("src/40_manuscript/main_figures/figure_1/figure_1f.R")

# Define input files
chronos_input_file <- glue::glue("{get_config()$input_dir}/name=chronos-pd_phase_1_odds_ratio_results/version=none/dataset.parquet")
us_input_file <- glue::glue("{get_config()$input_dir}/name=us_population_high_confidence_odds_ratio_results/version=none/dataset.parquet")

# Define output files
output_figure_file <- glue::glue("{get_config()$output_dir}/main_figures/figure_1f.pdf")
dir.create(dirname(output_figure_file), showWarnings = FALSE, recursive = TRUE)

# Load and combine data
data_list <- load_and_combine_odds_ratio_data(
  chronos_input_file |> arrow::open_dataset() |> dplyr::collect(),
  us_input_file |> arrow::open_dataset() |> dplyr::collect()
)

# Identify significant codes for pre-onset stage
significant_before_hc <- identify_significant_codes_for_stage(
  data_list$chronos_res,
  data_list$us_res,
  stage = "before"
)

# Prepare plotting data
plotting_data <- prepare_forest_plot_data(
  data_list$all_res,
  significant_before_hc,
  stage = "before"
)

# Create pre-onset plot
p <- create_forest_plot(
  plotting_data = plotting_data,
  title = "Pre-onset Comorbidities"
)

# Save plot
ggplot2::ggsave(
  filename = output_figure_file,
  plot = p,
  width = 13,
  height = 8,
  units = "in",
  dpi = 600,
  bg = "white"
)
