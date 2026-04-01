## SPDX-FileCopyrightText: 2026 Alkahest, Inc. <chronospd_data_access@alkahest.com>
##
## SPDX-License-Identifier: MIT

here::i_am("src/40_manuscript/main_figures/figure_2/figure_2h.R")

input_dataset_files <- list(
  "SomaScan" = glue::glue("{get_config()$input_dir}/name=chronos_somalogic_merged-pd_master_summary_2/version=29828/dataset.parquet"),
  "Explor HT" = glue::glue("{get_config()$input_dir}/name=chronos_olink_merged-pd_master_summary_2/version=29829/dataset.parquet"),
  "NULISA CNS" = glue::glue("{get_config()$input_dir}/name=chronos_alamar_merged-pd_master_summary_2/version=29830/dataset.parquet"),
  "Discovery MS" = glue::glue("{get_config()$input_dir}/name=chronos_biognosys_merged-pd_master_summary_2/version=29831/dataset.parquet")
)

output_figure_file <- glue::glue("{get_config()$output_dir}/main_figures/figure_2h.pdf")
dir.create(dirname(output_figure_file), showWarnings = FALSE, recursive = TRUE)

internal_protein_list_compare <- purrr::map(
  input_dataset_files,
  \(input_dataset_file) {
    arrow::open_dataset(input_dataset_file) |>
      dplyr::filter(q_val_all < 0.05) |>
      dplyr::distinct(uniprot_id) |>
      dplyr::pull(as_vector = TRUE)
  }
)

p <- VennDiagram::venn.diagram(
  x = internal_protein_list_compare,
  filename = NULL,
  output = TRUE,
  col = "transparent",
  fill = c("#00796B80", "#00B0F080", "#FFB30080", "#A2002580"), # , "#CCCCCC"
  cex = 1,
  fontface = "bold",
  cat.cex = 1,
  cat.fontface = "bold",
  margin = 0.05
)

ggplot2::ggsave(
  filename = output_figure_file,
  plot = p,
  width = 4.5, height = 4.5, units = "in",
  dpi = 600,
  bg = "white"
)
