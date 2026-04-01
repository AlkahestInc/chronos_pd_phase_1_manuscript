## SPDX-FileCopyrightText: 2026 Alkahest, Inc. <chronospd_data_access@alkahest.com>
##
## SPDX-License-Identifier: MIT

here::i_am("src/40_manuscript/main_figures/figure_2/figure_2i.R")

input_dataset_files <- list(
  "SomaScan" = glue::glue("{get_config()$input_dir}/name=chronos_somalogic_merged-pd_master_summary_2/version=29828/dataset.parquet"),
  "Explor HT" = glue::glue("{get_config()$input_dir}/name=chronos_olink_merged-pd_master_summary_2/version=29829/dataset.parquet"),
  "NULISA CNS" = glue::glue("{get_config()$input_dir}/name=chronos_alamar_merged-pd_master_summary_2/version=29830/dataset.parquet"),
  "Discovery MS" = glue::glue("{get_config()$input_dir}/name=chronos_biognosys_merged-pd_master_summary_2/version=29831/dataset.parquet")
)

paper_PD_list_freq <- arrow::open_dataset(
  glue::glue("{get_config()$input_dir}/name=paper_PD_list_freq/version=Nov2025/dataset.parquet")
) |>
  dplyr::collect()

output_figure_file <- glue::glue("{get_config()$output_dir}/main_figures/figure_2i.pdf")
dir.create(dirname(output_figure_file), showWarnings = FALSE, recursive = TRUE)

external_protein_list_compare <- purrr::map(
  input_dataset_files,
  \(input_dataset_file) {
    arrow::open_dataset(input_dataset_file) |>
      dplyr::filter(q_val_all < 0.05) |>
      dplyr::collect() |>
      tidyr::separate_longer_delim(
        uniprot_id,
        delim = stringr::regex("[\\|;:]")
      ) |>
      dplyr::distinct(uniprot_id) |>
      dplyr::pull()
  }
)
external_protein_list_compare[["Literature PD"]] <- paper_PD_list_freq |>
  dplyr::distinct(uniprot_id) |>
  dplyr::pull()

upset_data <- external_protein_list_compare |>
  UpSetR::fromList()

p <- upset_data |>
  UpSetR::upset(
    sets = c("SomaScan", "Explor HT", "NULISA CNS", "Discovery MS", "Literature PD"),
    order.by = "freq",
    text.scale = 2
  )

pdf(output_figure_file, width = 12, height = 6)
print(p)
dev.off()
