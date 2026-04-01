## SPDX-FileCopyrightText: 2026 Alkahest, Inc. <chronospd_data_access@alkahest.com>
##
## SPDX-License-Identifier: MIT

here::i_am("src/40_manuscript/supplementary_tables/table_s8.R")

input_dataset_files <- list(
  "Somalogic-11k" = glue::glue("{get_config()$input_dir}/name=chronos_somalogic_joint_model-pd_master_summary/version=29819/dataset.parquet"),
  "Olink-HT-5k" = glue::glue("{get_config()$input_dir}/name=chronos_olink_joint_model-pd_master_summary/version=29820/dataset.parquet"),
  "AlamarBio-CNS-120" = glue::glue("{get_config()$input_dir}/name=chronos_alamar_joint_model-pd_master_summary/version=29821/dataset.parquet")
)

# Define output file
output_table_file <- glue::glue("{get_config()$output_dir}/supplementary_tables/table_s8.tsv")
dir.create(dirname(output_table_file), showWarnings = FALSE, recursive = TRUE)


tbl <- purrr::map(
  input_dataset_files,
  ~ {
    .x |>
      arrow::open_dataset() |>
      dplyr::select(
        "Assay_ID" = "probe_id",
        "uniprot_id",
        "gene_symbol",
        "Estimate" = "estimate",
        "Standard Error" = "std_error",
        "lower_95",
        "upper_95",
        "Pvalue" = "p_val",
        "rhat",
        "Qvalue" = "q_val",
        "hr",
        "hr_lower",
        "hr_upper",
        "risk",
        "risk_lower",
        "risk_upper"
      ) |>
      dplyr::collect()
  }
) |>
  dplyr::bind_rows(.id = "platform") |>
  dplyr::arrange(
    dplyr::desc(abs(Estimate))
  )

tbl |>
  readr::write_tsv(output_table_file)
