## SPDX-FileCopyrightText: 2026 Alkahest, Inc. <chronospd_data_access@alkahest.com>
##
## SPDX-License-Identifier: MIT

here::i_am("src/40_manuscript/supplementary_tables/table_s4.R")

# Define output files
output_file <- glue::glue("{get_config()$output_dir}/supplementary_tables/table_s4.tsv")
dir.create(dirname(output_file), showWarnings = FALSE, recursive = TRUE)

spec <- tibble::tribble(
  ~file,                                                                                                                      ~proteomics_platform,     ~analysis_window,
  glue::glue("{get_config()$input_dir}/name=chronos_somalogic_all-pd_master_summary/version=29806/dataset.parquet"),          "Somalogic-11k",          "Pre- and post-onset",
  glue::glue("{get_config()$input_dir}/name=chronos_somalogic_pre_onset-pd_master_summary/version=29807/dataset.parquet"),    "Somalogic-11k",          "Pre-onset",
  glue::glue("{get_config()$input_dir}/name=chronos_somalogic_post_onset-pd_master_summary/version=29808/dataset.parquet"),   "Somalogic-11k",          "Post-onset",
  glue::glue("{get_config()$input_dir}/name=chronos_olink_all-pd_master_summary/version=29809/dataset.parquet"),              "Olink-HT-5k",            "Pre- and post-onset",
  glue::glue("{get_config()$input_dir}/name=chronos_olink_all_pre_onset-pd_master_summary/version=29810/dataset.parquet"),    "Olink-HT-5k",            "Pre-onset",
  glue::glue("{get_config()$input_dir}/name=chronos_olink_all_post_onset-pd_master_summary/version=29811/dataset.parquet"),   "Olink-HT-5k",            "Post-onset",
  glue::glue("{get_config()$input_dir}/name=chronos_alamar_all-pd_master_summary/version=29812/dataset.parquet"),             "AlamarBio-CNS-120",      "Pre- and post-onset",
  glue::glue("{get_config()$input_dir}/name=chronos_alamar_pre_onset-pd_master_summary/version=29813/dataset.parquet"),       "AlamarBio-CNS-120",      "Pre-onset",
  glue::glue("{get_config()$input_dir}/name=chronos_alamar_post_onset-pd_master_summary/version=29814/dataset.parquet"),      "AlamarBio-CNS-120",      "Post-onset",
  glue::glue("{get_config()$input_dir}/name=chronos_biognosys_all-pd_master_summary/version=29864/dataset.parquet"),          "Discovery MS",           "Pre- and post-onset",
  glue::glue("{get_config()$input_dir}/name=chronos_biognosys_pre_onset-pd_master_summary/version=29866/dataset.parquet"),    "Discovery MS",           "Pre-onset",
  glue::glue("{get_config()$input_dir}/name=chronos_biognosys_post_onset-pd_master_summary/version=29867/dataset.parquet"),   "Discovery MS",           "Post-onset",
)

process_model_group_results <- function(file, proteomics_platform, analysis_window) {
  arrow::open_dataset(file) |>
    dplyr::collect() |>
    dplyr::mutate(
      "Assay ID" = probe_id,
      "UniProt ID" = uniprot_id,
      "Gene Symbol" = gene_symbol,
      "Estimate" = estimate,
      "Standard Error" = std_err,
      "Cohen's d" = cohens_d,
      "Pvalue" = p_val,
      "Qvalue" = q_val,
      "Proteomics Platform" = proteomics_platform,
      "Analysis window" = analysis_window,
      .keep = "none"
    )
}

result <- purrr::pmap_dfr(spec, process_model_group_results) |>
  dplyr::arrange(Pvalue)

result |>
  readr::write_tsv(output_file)
