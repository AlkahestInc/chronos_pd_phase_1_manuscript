## SPDX-FileCopyrightText: 2026 Alkahest, Inc. <chronospd_data_access@alkahest.com>
##
## SPDX-License-Identifier: MIT

here::i_am("src/40_manuscript/supplementary_tables/table_s10.R")

# Input paths
dir_path <- glue::glue("{get_config()$checkpoint_dir}/joint_modeling/comorbidity/other_comorbidity/summary_only")
OR_phase1_chronic_ref_path <- glue::glue("{get_config()$checkpoint_dir}/joint_modeling/comorbidity/significant_chronic_ref.parquet")

# Output paths
output_file <- glue::glue("{get_config()$output_dir}/supplementary_tables/table_s10.tsv")

# ==============================================================================
# Load and process SomaLogic data
# ==============================================================================

files <- list.files(
  dir_path,
  pattern    = "^somalogic_JM\\..*\\.parquet$",
  full.names = TRUE
)

combined_soma <- tibble::tibble(file = files) |>
  dplyr::mutate(file_name = basename(file)) |>
  tidyr::extract(
    col   = file_name,
    into  = c("platform", "diagnosis_code_3", "protein_index", "protein_from_name"),
    regex = "^(somalogic_JM)\\.([^\\.]+)\\.([0-9]+)\\.(.*)\\.parquet$"
  ) |>
  dplyr::mutate(
    summary_tbl = purrr::map(file, function(x){
      x |>
        arrow::open_dataset() |>
        dplyr::collect()
    })
  ) |>
  dplyr::select(-protein_from_name) |>
  tidyr::unnest(summary_tbl) |>
  dplyr::select(-file) |>
  dplyr::mutate(
    combination = paste0(diagnosis_code_3, "_", protein)
  )

dim(combined_soma) # 831  15

# ==============================================================================
# Load and process Olink data
# ==============================================================================

files <- list.files(
  dir_path,
  pattern    = "^olink_JM\\..*\\.parquet$",
  full.names = TRUE
)

combined_olink <- tibble::tibble(file = files) |>
  dplyr::mutate(file_name = basename(file)) |>
  tidyr::extract(
    col   = file_name,
    into  = c("platform", "diagnosis_code_3", "protein_index", "protein_from_name"),
    regex = "^(olink_JM)\\.([^\\.]+)\\.([0-9]+)\\.(.*)\\.parquet$"
  ) |>
  dplyr::mutate(
    summary_tbl = purrr::map(file, function(x){
      x |>
        arrow::open_dataset() |>
        dplyr::collect()
    })
  ) |>
  dplyr::select(-protein_from_name) |>
  tidyr::unnest(summary_tbl) |>
  dplyr::select(-file) |>
  dplyr::mutate(
    combination = paste0(diagnosis_code_3, "_", protein)
  )

dim(combined_olink) # 571  15; one missing

# ==============================================================================
# Load and process Alamar data
# ==============================================================================

files <- list.files(
  dir_path,
  pattern    = "^alamar_JM\\..*\\.parquet$",
  full.names = TRUE
)

combined_alamar <- tibble::tibble(file = files) |>
  dplyr::mutate(file_name = basename(file)) |>
  tidyr::extract(
    col   = file_name,
    into  = c("platform", "diagnosis_code_3", "protein_index", "protein_from_name"),
    regex = "^(alamar_JM)\\.([^\\.]+)\\.([0-9]+)\\.(.*)\\.parquet$"
  ) |>
  dplyr::mutate(
    summary_tbl = purrr::map(file, function(x){
      x |>
        arrow::open_dataset() |>
        dplyr::collect()
    })
  ) |>
  dplyr::select(-protein_from_name) |>
  tidyr::unnest(summary_tbl) |>
  dplyr::select(-file) |>
  dplyr::mutate(
    combination = paste0(diagnosis_code_3, "_", protein)
  )

dim(combined_alamar) # 468  15

# ==============================================================================
# Load reference data
# ==============================================================================

OR_phase1_chronic_ref <- OR_phase1_chronic_ref_path |>
  arrow::open_dataset() |>
  dplyr::collect()

# ==============================================================================
# Combine all platforms and join with reference data
# ==============================================================================

JM_condition_results <- combined_soma |>
  dplyr::bind_rows(combined_olink) |>
  dplyr::bind_rows(combined_alamar) |>
  dplyr::left_join(
    protein_list |>
      dplyr::select(-platform),
    by = c("protein" = "probe_id")
  ) |>
  dplyr::left_join(
    OR_phase1_chronic_ref,
    by = "diagnosis_code_3"
  ) |>
  dplyr::relocate(ICD10_description)

dim(JM_condition_results) # 1870   19

#saveRDS(JM_condition_results, "/workbench/ying_folder/Work_for_Manuscript/phase1/JM/comorbidity_ICDconverted/nonPD_conditions/JM_condition_results.RDS")

# ==============================================================================
# Filter conditions and calculate hazard ratios and risk metrics
# ==============================================================================

JM_condition_results_saved <- JM_condition_results |>
  dplyr::filter(ICD10_description %nin% c("Essential (primary) hypertension", "Type 2 diabetes mellitus")) |>
  dplyr::mutate(
    # Calculate hazard ratios
    hr       = exp(estimate),
    hr_lower = exp(lower_95),
    hr_upper = exp(upper_95),
    # Calculate risk per 0.1 unit change
    risk       = ifelse(estimate > 0, exp(estimate * 0.1) - 1, 1 - exp(estimate * 0.1)),
    risk_lower = ifelse(estimate > 0, exp(lower_95 * 0.1) - 1, 1 - exp(lower_95 * 0.1)),
    risk_upper = ifelse(estimate > 0, exp(upper_95 * 0.1) - 1, 1 - exp(upper_95 * 0.1))
  ) |>
  dplyr::select(
    protein, uniprot_id, gene_symbol, diagnosis_code_3, ICD10_description,
    term, estimate, std_error, lower_95, upper_95, p_value, rhat, dic, waic, lpml,
    hr, hr_lower, hr_upper, risk, risk_lower, risk_upper, platform
  ) |>
  dplyr::mutate(
    platform = dplyr::case_when(
      platform == "somalogic_JM" ~ "SomaScanV5",
      platform == "olink_JM"     ~ "OlinkHT",
      TRUE                       ~ "Alamarbio_CNS"
    )
  ) |>
  dplyr::rename(
    probe_id = protein,
    p_val    = p_value
  ) |>
  # Select and rename columns for publication
  dplyr::select(
    `Assay ID` = probe_id,
    `UniProt ID` = uniprot_id,
    `Gene Symbol` = gene_symbol,
    `Platform` = platform,
    `ICD-10 Code 3 digit` = diagnosis_code_3,
    `ICD-10 Description` = ICD10_description,
    `Estimate` = estimate,
    `Standard Error` = std_error,
    `Hazard Ratio` = hr,
    `HR Lower 95% CI` = hr_lower,
    `HR Upper 95% CI` = hr_upper,
    `P-value` = p_val,
    `Rhat` = rhat,
    `DIC` = dic,
    `WAIC` = waic,
    `LPML` = lpml,
    `Risk % per 0.1 SD` = risk,
    `Risk Lower 95% CI` = risk_lower,
    `Risk Upper 95% CI` = risk_upper
  ) |>
  # Clean up platform names
  dplyr::mutate(
    Platform = Platform |> dplyr::recode_values(
      "SomaScanV5" ~ "Somalogic-11k",
      "OlinkHT" ~ "Olink-HT-5k",
      "Alamarbio_CNS" ~ "AlamarBio-CNS-120"
    )
  ) |>
  dplyr::arrange(`P-value`)

dim(JM_condition_results_saved) # 1798   22

# ==============================================================================
# Save output table
# ==============================================================================

JM_condition_results_saved |>
  readr::write_tsv(output_file)
