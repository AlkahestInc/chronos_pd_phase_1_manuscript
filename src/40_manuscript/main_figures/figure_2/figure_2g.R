## SPDX-FileCopyrightText: 2026 Alkahest, Inc. <chronospd_data_access@alkahest.com>
##
## SPDX-License-Identifier: MIT

here::i_am("src/40_manuscript/main_figures/figure_2/figure_2g.R")

input_dataset_files <- list(
  "SomaScan" = glue::glue("{get_config()$input_dir}/name=chronos_somalogic_merged-pd_master_summary_2/version=29828/dataset.parquet"),
  "Explor HT" = glue::glue("{get_config()$input_dir}/name=chronos_olink_merged-pd_master_summary_2/version=29829/dataset.parquet"),
  "NULISA CNS" = glue::glue("{get_config()$input_dir}/name=chronos_alamar_merged-pd_master_summary_2/version=29830/dataset.parquet"),
  "Discovery MS" = glue::glue("{get_config()$input_dir}/name=chronos_biognosys_merged-pd_master_summary_2/version=29831/dataset.parquet")
)

output_figure_file <- glue::glue("{get_config()$output_dir}/main_figures/figure_2g.pdf")
dir.create(dirname(output_figure_file), showWarnings = FALSE, recursive = TRUE)

external_protein_list_compare <- purrr::map(
  input_dataset_files,
  \(input_dataset_file) {
    arrow::open_dataset(input_dataset_file) |>
      dplyr::distinct(uniprot_id, gene_symbol, p_val_all, log2_fc_all, q_val_all) |>
      dplyr::collect()
  }
)

heatmap_data <- external_protein_list_compare |>
  dplyr::bind_rows(.id = "platform") |>
  dplyr::mutate(
    signed = -log10(p_val_all) * sign(log2_fc_all)
  ) |>
  #####
  # Where a uniprot/gene_symbol has more than 1 measurement on a platform,
  # keep only the best (according to p_val_all) one
  #####
  dplyr::group_by(uniprot_id, gene_symbol, platform) |>
  dplyr::slice_min(
    p_val_all,
    n = 1,
    with_ties = TRUE
  ) |>
  dplyr::ungroup() |>
  #####
  # Remove uniprot_id's that don't have any q_val_all < 0.05
  ####
  dplyr::group_by(uniprot_id) |>
  dplyr::filter(
    min(q_val_all) < 0.05
  ) |>
  dplyr::ungroup() |>
  #####
  # Reshape the data
  #####
  dplyr::select(
    uniprot_id,
    gene_symbol,
    platform,
    signed
  ) |>
  tidyr::pivot_wider(
    names_from  = platform,
    values_from = signed
  ) |>
  dplyr::mutate(
    id = paste(uniprot_id, gene_symbol, sep = "_")
  ) |>
  dplyr::select(
    id,
    where(is.double)
  )

# Generate a dendrogram for row ordering of the heatmap
row_dend <- heatmap_data |>
  dplyr::mutate(
    dplyr::across(tidyselect::everything(), ~ tidyr::replace_na(.x, 0))
  ) |>
  tibble::column_to_rownames("id") |>
  as.matrix() |>
  dist() |>
  hclust()

p <- heatmap_data |>
  tibble::column_to_rownames("id") |>
  as.matrix() |>
  ComplexHeatmap::Heatmap(
    cluster_rows = row_dend,
    column_order = c("SomaScan", "Explor HT", "Discovery MS", "NULISA CNS"),
    show_row_names = FALSE,
    name = "signed\n-log10(p)",
    na_col = "grey80",
    column_names_rot = 45,
    col = circlize::colorRamp2(
      c(-9, 0, 11),
      c("blue", "white", "red")
    )
  )

pdf(file = output_figure_file, width = 3.0, height = 5.2)
ComplexHeatmap::draw(p)
dev.off()
