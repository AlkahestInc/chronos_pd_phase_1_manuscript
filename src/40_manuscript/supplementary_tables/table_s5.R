## SPDX-FileCopyrightText: 2026 Alkahest, Inc. <chronospd_data_access@alkahest.com>
##
## SPDX-License-Identifier: MIT

here::i_am("src/40_manuscript/supplementary_tables/table_s5.R")

input_dataset_files <- list(
  "Somalogic-11k" = glue::glue("{get_config()$input_dir}/name=chronos_somalogic_all-pd_master_summary/version=29806/dataset.parquet"),
  "Olink-HT-5k" = glue::glue("{get_config()$input_dir}/name=chronos_olink_all-pd_master_summary/version=29809/dataset.parquet"),
  "AlamarBio-CNS-120" = glue::glue("{get_config()$input_dir}/name=chronos_alamar_all-pd_master_summary/version=29812/dataset.parquet")
)

# Define output file
output_table_file <- glue::glue("{get_config()$output_dir}/supplementary_tables/table_s5.tsv")
dir.create(dirname(output_table_file), showWarnings = FALSE, recursive = TRUE)

PD_proteins <- purrr::map(
  input_dataset_files,
  ~ {
    .x |>
      arrow::open_dataset() |>
      dplyr::filter(
        q_val < 0.05
      ) |>
      dplyr::select(gene_symbol) |>
      dplyr::collect()
  }
)

pairs <- combn(names(PD_proteins), 2, simplify = FALSE)

pairwise_intersections <- purrr::map(
  combn(names(PD_proteins), 2, simplify = FALSE),
  ~ {
    tibble::tibble(
      platform1 = .x[1],
      platform2 = .x[2],
      Intersect = glue::glue("{platform1} U {platform2}"),
      gene_symbol = intersect(
        PD_proteins[[.x[1]]]$gene_symbol,
        PD_proteins[[.x[2]]]$gene_symbol
      )
    )
  }
) |>
  dplyr::bind_rows()

tbl <- pairwise_intersections |>
  dplyr::group_by(
    Intersect
  ) |>
  dplyr::summarise(
    "Discovery MS" = NA,
    "Number of proteins" = dplyr::n_distinct(gene_symbol),
    "Proteins" = paste(sort(unique(gene_symbol)), collapse = ", "),
    .groups = "drop"
  ) |>
  dplyr::arrange(
    dplyr::desc(`Number of proteins`)
  )

tbl |>
  readr::write_tsv(output_table_file)
