## SPDX-FileCopyrightText: 2026 Alkahest, Inc. <chronospd_data_access@alkahest.com>
##
## SPDX-License-Identifier: MIT

here::i_am("src/20_analyses/non_linear_mixed/alamar.R")

#####
# Define inputs
#####
alamar_tibble_file <- glue::glue("{get_config()$input_dir}/name=Chronos_phase1_proteomics_alamar_analysis_ready_three_tibble/version=29649/dataset.rds")
alamar_all_file <- glue::glue("{get_config()$input_dir}/name=chronos_alamar_all-pd_master_summary/version=29812/dataset.parquet")
alamar_model_file <- glue::glue("{get_config()$input_dir}/name=alamar_screen_model_comparison/version=none/dataset.rds")

#####
# Define outputs
#####
outfile_ls <- glue::glue("{get_config()$checkpoint_dir}/non_linear_mixed/alamar/alamar_all_ls.rds")
outfile_nc <- glue::glue("{get_config()$checkpoint_dir}/non_linear_mixed/alamar/alamar_all_nc.rds")
dir.create(dirname(outfile_ls), showWarnings = FALSE, recursive = TRUE)

#####
# Analyses
#####
alamar_tibble <- alamar_tibble_file |>
  readr::read_rds()

alamar_table <- alamar_tibble$meta_data |>
  dplyr::inner_join(
    alamar_tibble$prot_data,
    by = "ID"
  )

alamar_all <- alamar_all_file |>
  arrow::open_dataset() |>
  dplyr::collect()
alamar_all_hits <- unique(alamar_all$probe_id[which(alamar_all$q_val < 0.05)]) # 14
alamar_all_hits

# alamar selected model structure:
alamar_model <- alamar_model_file |>
  readr::read_rds()

alamar_all_ls <- run_linear_spline_models_pb(
  data = alamar_table,
  proteins = alamar_all_hits,
  lme_results = alamar_model,
  knot_range = -7:5,
  mc.cores = parallelly::availableCores()
)

readr::write_rds(alamar_all_ls, outfile_ls)

selected_proteins_BIC <- names(alamar_all_ls)[
  sapply(alamar_all_ls, function(df) any(df$diff_lsBIC > 2, na.rm = TRUE))
]

if (length(selected_proteins_BIC) > 0) {
  alamar_all_ls_selected_BIC <- dplyr::bind_rows(alamar_all_ls[selected_proteins_BIC], .id = "protein") |>
    dplyr::filter(diff_lsBIC > 2 & group_diff_sig != "neither")

  alamar_all_ls_selected_BIC |>
    dplyr::group_by(protein) |>
    dplyr::filter(diff_lsBIC == max(diff_lsBIC))

  alamar_all_ls_selected_BIC |>
    dplyr::filter(group_diff_sig != "neither") |>
    dplyr::distinct(protein)
} else {
  cat("No alamar proteins with diff_lsBIC > 2\n")
}

alamar_all_nc <- run_ncspline_models_pb(
  data = alamar_table,
  proteins = alamar_all_hits,
  lme_results = alamar_model,
  df_range = 2:6,
  mc.cores = parallelly::availableCores()
)

readr::write_rds(alamar_all_nc, outfile_nc)

selected_proteins_BIC <- names(alamar_all_nc)[
  sapply(alamar_all_nc, function(df) any(df$diff_lmeBIC > 2, na.rm = TRUE))
]

if (length(selected_proteins_BIC) > 0) {
  dplyr::bind_rows(alamar_all_nc[selected_proteins_BIC], .id = "protein")
} else {
  cat("No alamar proteins with diff_lmeBIC > 2 for NC spline\n")
}
