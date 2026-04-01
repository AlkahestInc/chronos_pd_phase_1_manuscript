## SPDX-FileCopyrightText: 2026 Alkahest, Inc. <chronospd_data_access@alkahest.com>
##
## SPDX-License-Identifier: MIT

here::i_am("src/20_analyses/non_linear_mixed/somalogic.R")

#####
# Define inputs
#####
somalogic_tibble_file <- glue::glue("{get_config()$input_dir}/name=Chronos_phase1_proteomics_somalogic_analysis_ready_three_tibble/version=29650/dataset.rds")
somalogic_all_file <- glue::glue("{get_config()$input_dir}/name=chronos_somalogic_all-pd_master_summary/version=29806/dataset.parquet")
somalogic_model_file <- glue::glue("{get_config()$input_dir}/name=somalogic_screen_model_comparison/version=none/dataset.rds")

#####
# Define outputs
#####
outfile_ls <- glue::glue("{get_config()$checkpoint_dir}/non_linear_mixed/somalogic/somalogic_all_ls.rds")
outfile_nc <- glue::glue("{get_config()$checkpoint_dir}/non_linear_mixed/somalogic/somalogic_all_nc.rds")
dir.create(dirname(outfile_ls), showWarnings = FALSE, recursive = TRUE)

#####
# Analyses
#####
somalogic_tibble <- somalogic_tibble_file |>
  readr::read_rds()

somalogic_table <- somalogic_tibble$meta_data |>
  dplyr::inner_join(
    somalogic_tibble$prot_data,
    by = "ID"
  )

somalogic_all <- somalogic_all_file |>
  arrow::open_dataset() |>
  dplyr::collect()
somalogic_all_hits <- unique(somalogic_all$probe_id[which(somalogic_all$q_val < 0.05)]) # 417
somalogic_all_hits

# somalogic selected model structure:
somalogic_model <- somalogic_model_file |>
  readr::read_rds()

somalogic_all_ls <- run_linear_spline_models_pb(
  data = somalogic_table,
  proteins = somalogic_all_hits,
  lme_results = somalogic_model,
  knot_range = -7:5,
  mc.cores = parallelly::availableCores()
)

readr::write_rds(somalogic_all_ls, outfile_ls)

selected_proteins_BIC <- names(somalogic_all_ls)[
  sapply(somalogic_all_ls, function(df) any(df$diff_lsBIC > 2, na.rm = TRUE))
]

if (length(selected_proteins_BIC) > 0) {
  somalogic_all_ls_selected_BIC <- dplyr::bind_rows(somalogic_all_ls[selected_proteins_BIC], .id = "protein") |>
    dplyr::filter(diff_lsBIC > 2 & group_diff_sig != "neither")

  somalogic_all_ls_selected_BIC |>
    dplyr::group_by(protein) |>
    dplyr::filter(diff_lsBIC == max(diff_lsBIC))

  somalogic_all_ls_selected_BIC |>
    dplyr::filter(group_diff_sig != "neither") |>
    dplyr::distinct(protein)
} else {
  cat("No somalogic proteins with diff_lsBIC > 2\n")
}

somalogic_all_nc <- run_ncspline_models_pb(
  data = somalogic_table,
  proteins = somalogic_all_hits,
  lme_results = somalogic_model,
  df_range = 2:6,
  mc.cores = parallelly::availableCores()
)

readr::write_rds(somalogic_all_nc, outfile_nc)

selected_proteins_BIC <- names(somalogic_all_nc)[
  sapply(somalogic_all_nc, function(df) any(df$diff_lmeBIC > 2, na.rm = TRUE))
]

if (length(selected_proteins_BIC) > 0) {
  dplyr::bind_rows(somalogic_all_nc[selected_proteins_BIC], .id = "protein")
} else {
  cat("No somalogic proteins with diff_lmeBIC > 2 for NC spline\n")
}
