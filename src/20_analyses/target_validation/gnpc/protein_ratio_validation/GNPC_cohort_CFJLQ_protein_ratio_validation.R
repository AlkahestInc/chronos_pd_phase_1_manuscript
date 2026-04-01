## SPDX-FileCopyrightText: 2026 Alkahest, Inc. <chronospd_data_access@alkahest.com>
##
## SPDX-License-Identifier: MIT

here::i_am("src/20_analyses/target_validation/gnpc/protein_ratio_validation/GNPC_cohort_CFJLQ_protein_ratio_validation.R")

input_dir <- "/home/rstudio"

# load GNPC annotation tibble
annotations <- glue::glue("{input_dir}/data_prep/annotation/annotation_tbl_v6.rds") |>
  readr::read_rds()

# load protein pair biomarkers
prot_ratio_tibble <- glue::glue(
  "{get_config()$input_dir}/name={name}/version={version}/dataset.parquet",
  name = "Chronos_phase1_protein_ratio_list_1",
  version = 29452
) |>
  arrow::open_dataset() |>
  dplyr::collect()

# probe names from somalogic starts with `seq`, from olink starts with `OID`, and alamarbio directly given gene symbols
protein_pair_tibble_list <- list(
  chronos_PD_ratio_list1 = prot_ratio_tibble
)

# check/create output path
output_dir <- glue::glue("{get_config()$checkpoint_dir}/GNPC_validation/model/protein_ratio_validation/protein_ratio_validation_result")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Build "other" analyte mapping only if missing
# other_analyte_mapping Data frame for other platform mappings.
mapping_tibbles <- build_mapping_tibbles(protein_pair_tibble_list)
other_analyte_mapping <- mapping_tibbles |>
  purrr::reduce(dplyr::bind_rows) |>
  dplyr::distinct(.keep_all = TRUE)

for (selected_contributor_code in c("C", "F", "J_plasma", "L", "Q_plasma")) {
  cat("For GNPC cohort : ", selected_contributor_code, "\n")

  if (selected_contributor_code == "J_plasma") {
    matching_ratio_str <- "1to1"
  } else {
    matching_ratio_str <- "2to1"
  }
  if (selected_contributor_code %in% c("C", "L", "Q_plasma")) {
    visit_str <- ""
  } else {
    visit_str <- "_last_visit" # "_first_visit"
  }

  # load gnpc protein intensity and other data
  metadata <- glue::glue("{get_config()$checkpoint_dir}/data_prep/gnpc_cohort_cleanup_matched_data/cohort_{selected_contributor_code}{visit_str}_outliers_rm_{matching_ratio_str}_matched_grp_filtered_metadata.rds") |>
    readr::read_rds()
  protein_intensity_log2 <- glue::glue("{get_config()$checkpoint_dir}/data_prep/gnpc_cohort_cleanup_matched_data/cohort_{selected_contributor_code}{visit_str}_outliers_rm_{matching_ratio_str}_matched_grp_filtered_protein_intensity_log2.rds") |>
    readr::read_rds()

  all_results <- list()
  all_summary_stats <- list()
  all_significant_ratios <- list()
  for (tbl_name in names(protein_pair_tibble_list)) {
    cat("Calculating for protein pair table : ", tbl_name, "\n")
    results <- analyze_multiple_ratio_by_tibble(
      protein_pair_tibble = protein_pair_tibble_list[[tbl_name]],
      protein_ratio_metadata = metadata,
      protein_intensity_log2 = protein_intensity_log2
    ) |>
      dplyr::mutate(
        protein_pair_tbl_name = !!tbl_name
      ) |>
      dplyr::relocate(protein_pair_tbl_name)

    # adjusted logistic regression for single protein ratio
    results_summary <- summarize_ratio_results(results)

    # Collect
    all_results[[tbl_name]] <- results

    summary_stats <- results_summary$summary_stats |>
      dplyr::mutate(
        protein_pair_tbl_name = !!tbl_name
      ) |>
      dplyr::relocate(protein_pair_tbl_name)
    all_summary_stats[[tbl_name]] <- summary_stats

    significant_ratios <- results_summary$significant_ratios |>
      dplyr::mutate(
        protein_pair_tbl_name = !!tbl_name
      ) |>
      dplyr::relocate(protein_pair_tbl_name)
    all_significant_ratios[[tbl_name]] <- significant_ratios

    cat("Done for for protein pair table : ", tbl_name, " ; results dims (", dim(results), ")\n")
  }

  output_file_path <- glue::glue("{output_dir}/cohort_{selected_contributor_code}_outliers_rm_{matching_ratio_str}_matched_grp_3tbl_protein_ratio_validation_results.csv")
  output_summary_stats_path <- glue::glue("{output_dir}/cohort_{selected_contributor_code}_outliers_rm_{matching_ratio_str}_matched_grp_3tbl_protein_ratio_validation_result_summary_stats.csv")
  output_summary_significant_ratios_path <- glue::glue("{output_dir}/cohort_{selected_contributor_code}_outliers_rm_{matching_ratio_str}_matched_grp_3tbl_protein_ratio_validation_result_summary_significant_ratios.csv")

  final_results <- do.call(rbind, all_results)
  # filter out unmatched or unconverged rows (NA rows)
  final_results <- final_results |> dplyr::filter(converged)

  final_results |>
    readr::write_csv(output_file_path)
  cat("Output all protein ratio validation result table for GNPC cohort : ", selected_contributor_code, "; results dims (", dim(final_results), ")\n")

  final_significant_ratios <- do.call(rbind, all_significant_ratios)
  final_significant_ratios |>
    readr::write_csv(output_summary_significant_ratios_path)
  cat("Output significant protein ratio table for GNPC cohort : ", selected_contributor_code, "; results dims (", dim(final_significant_ratios), ")\n")

  final_summary_stats <- do.call(rbind, all_summary_stats)
  final_summary_stats |>
    readr::write_csv(output_summary_stats_path)
  print(final_summary_stats)
  cat("Output protein ratio final summary table for GNPC cohort : ", selected_contributor_code, "\n")
  cat("Done for GNPC cohort : ", selected_contributor_code, "\n")
}
