## SPDX-FileCopyrightText: 2026 Alkahest, Inc. <chronospd_data_access@alkahest.com>
##
## SPDX-License-Identifier: MIT

here::i_am("src/20_analyses/non_linear_mixed/olink.R")

#####
# Define inputs
#####
olink_tibble_file <- glue::glue("{get_config()$input_dir}/name=Chronos_phase1_proteomics_olink_analysis_ready_three_tibble/version=29465/dataset.rds")
olink_all_file <- glue::glue("{get_config()$input_dir}/name=chronos_olink_all-pd_master_summary/version=29809/dataset.parquet")
olink_model_file <- glue::glue("{get_config()$input_dir}/name=olink_screen_model_comparison/version=none/dataset.rds")

#####
# Define outputs
#####
outfile_ls <- glue::glue("{get_config()$checkpoint_dir}/non_linear_mixed/olink/olink_all_ls.rds")
outfile_nc <- glue::glue("{get_config()$checkpoint_dir}/non_linear_mixed/olink/olink_all_nc.rds")
dir.create(dirname(outfile_ls), showWarnings = FALSE, recursive = TRUE)

#####
# Analyses
#####
olink_tibble <- olink_tibble_file |>
  readr::read_rds()

olink_table <- olink_tibble$meta_data |>
  dplyr::inner_join(
    olink_tibble$prot_data,
    by = "ID"
  )

olink_all <- olink_all_file |>
  arrow::open_dataset() |>
  dplyr::collect()
olink_all_hits <- unique(olink_all$probe_id[which(olink_all$q_val < 0.05)]) # 101
olink_all_hits

# olink selected model structure:
olink_model <- olink_model_file |>
  readr::read_rds()

olink_all_ls <- run_linear_spline_models_pb(
  data = olink_table,
  proteins = olink_all_hits,
  lme_results = olink_model,
  knot_range = -7:5,
  mc.cores = parallelly::availableCores()
)

readr::write_rds(olink_all_ls, outfile_ls)

selected_proteins_BIC <- names(olink_all_ls)[
  sapply(olink_all_ls, function(df) any(df$diff_lsBIC > 2, na.rm = TRUE))
]

if (length(selected_proteins_BIC) > 0) {
  olink_all_ls_selected_BIC <- dplyr::bind_rows(olink_all_ls[selected_proteins_BIC], .id = "protein") |>
    dplyr::filter(diff_lsBIC > 2 & group_diff_sig != "neither")

  olink_all_ls_selected_BIC |>
    dplyr::group_by(protein) |>
    dplyr::filter(diff_lsBIC == max(diff_lsBIC))

  olink_all_ls_selected_BIC |>
    dplyr::filter(group_diff_sig != "neither") |>
    dplyr::distinct(protein)
} else {
  cat("No olink proteins with diff_lsBIC > 2\n")
}

olink_all_nc <- run_ncspline_models_pb(
  data = olink_table,
  proteins = olink_all_hits,
  lme_results = olink_model,
  df_range = 2:6,
  mc.cores = parallelly::availableCores()
)

readr::write_rds(olink_all_nc, outfile_nc)

selected_proteins_BIC <- names(olink_all_nc)[
  sapply(olink_all_nc, function(df) any(df$diff_lmeBIC > 2, na.rm = TRUE))
]

if (length(selected_proteins_BIC) > 0) {
  dplyr::bind_rows(olink_all_nc[selected_proteins_BIC], .id = "protein")
} else {
  cat("No olink proteins with diff_lmeBIC > 2 for NC spline\n")
}
