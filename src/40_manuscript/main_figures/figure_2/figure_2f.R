## SPDX-FileCopyrightText: 2026 Alkahest, Inc. <chronospd_data_access@alkahest.com>
##
## SPDX-License-Identifier: MIT

here::i_am("src/40_manuscript/main_figures/figure_2/figure_2f.R")

input_dataset_file <- glue::glue("{get_config()$input_dir}/name=chronos_biognosys_merged-pd_master_summary_2/version=29831/dataset.parquet")
output_figure_file <- glue::glue("{get_config()$output_dir}/main_figures/figure_2f.pdf")
paper_pd_list_file <- glue::glue("{get_config()$input_dir}/name=paper_PD_list_freq/version=Nov2025/dataset.parquet")

# Create output directory if it doesn't exist
dir.create(dirname(output_figure_file), showWarnings = FALSE, recursive = TRUE)

# Load paper PD list
paper_pd_list <- arrow::open_dataset(paper_pd_list_file) |>
  dplyr::collect()

# Load and process dataset
dataset_result <- load_volcano_dataset(input_dataset_file, paper_pd_list)

# Create plot
p <- create_volcano_plot(dataset_result, paper_pd_list)

# Save plot
ggplot2::ggsave(
  filename = output_figure_file,
  plot = p,
  width = 4.3,
  height = 5.2,
  units = "in",
  dpi = 600,
  bg = "white"
)
