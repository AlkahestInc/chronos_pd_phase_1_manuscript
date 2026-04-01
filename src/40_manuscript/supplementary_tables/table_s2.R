## SPDX-FileCopyrightText: 2026 Alkahest, Inc. <chronospd_data_access@alkahest.com>
##
## SPDX-License-Identifier: MIT

here::i_am("src/40_manuscript/supplementary_tables/table_s2.R")

# Define inputs
chronos_input_file <- glue::glue("{get_config()$input_dir}/name=chronos-pd_phase_1_odds_ratio_results/version=none/dataset.parquet")
us_input_file <- glue::glue("{get_config()$input_dir}/name=us_population_high_confidence_odds_ratio_results/version=none/dataset.parquet")

# Define outputs
output_file_comorbiditiy_profiles <- glue::glue("{get_config()$output_dir}/supplementary_tables/table_s2/table_s2.tsv")
dir.create(dirname(output_file_comorbiditiy_profiles), showWarnings = FALSE, recursive = TRUE)

# Load data, select relevant columns, and filter for significance
chronos_sig_results <- chronos_input_file |>
  arrow::open_dataset() |>
  dplyr::select(
    `Analysis window` = stage,
    `ICD10 code (3 digits)` = DIAGNOSIS_CODE_CATEGORY,
    `ICD10 code description` = label,
    Estimate,
    `Standard Error` = StdError,
    Pvalue = PValue,
    Qvalue = FDR,
    OR,
    LowerCI,
    UpperCI
  ) |>
  dplyr::mutate(
    `Analysis window` = dplyr::case_when(
      `Analysis window` == "before" ~ "Pre-onset",
      `Analysis window` == "after" ~ "Post-onset"
    )
  ) |>
  dplyr::filter(
    Qvalue <= 0.05
  ) |>
  dplyr::arrange(
    dplyr::desc(`Analysis window`),
    Pvalue
  ) |>
  dplyr::collect()

# Save output
chronos_sig_results |>
  readr::write_tsv(output_file_comorbiditiy_profiles)
