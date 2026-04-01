## SPDX-FileCopyrightText: 2026 Alkahest, Inc. <chronospd_data_access@alkahest.com>
##
## SPDX-License-Identifier: MIT

here::i_am("src/20_analyses/target_validation/gnpc/single_protein_validation/GNPC_cohort_CFJLQ_single_protein_validation.R")

input_dir <- "/home/rstudio"

# load GNPC annotation tibble
annotations <- glue::glue("{input_dir}/data_prep/annotation/annotation_tbl_v6.rds") |>
  readr::read_rds()

# check/create output path
output_dir <- glue::glue("{get_config()$checkpoint_dir}/GNPC_validation/model/single_protein_validation/result")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# load chronos single protein biomarkers and filter by q value cutoff
q_val_cutoff <- 0.05
chronos_somalogic_protein_tibble <- glue::glue("{input_dir}/input_data/single_protein_validation/chronos_somalogic_protein_list_all.rds") |>
  readr::read_rds() |>
  dplyr::filter(
    q_val < q_val_cutoff
  )

chronos_olink_protein_tibble <- glue::glue("{input_dir}/input_data/single_protein_validation/chronos_olink_protein_list_all.rds") |>
  readr::read_rds() |>
  dplyr::filter(
    q_val < q_val_cutoff
  )

chronos_alamar_protein_tibble <- glue::glue("{input_dir}/input_data/single_protein_validation/chronos_alamar_protein_list_all.rds") |>
  readr::read_rds() |>
  dplyr::filter(
    q_val < q_val_cutoff
  )

# probe names from somalogic starts with `seq`, from olink starts with `OID`, and alamarbio directly given gene symbols
single_protein_tibble_list <- list(
  chronos_somalogic_protein_list = chronos_somalogic_protein_tibble,
  chronos_olink_protein_list = chronos_olink_protein_tibble,
  chronos_alamar_protein_list = chronos_alamar_protein_tibble
)

# Build "other" analyte mapping only if missing
# other_analyte_mapping Data frame for other platform mappings.
mapping_tibbles <- list(
  chronos_somalogic_protein_tibble |>
    dplyr::select(
      platform_id = platform_id,
      uniprot_id  = uniprot_id,
      gene_symbol = gene_symbol
    ),
  chronos_olink_protein_tibble |>
    dplyr::select(
      platform_id = platform_id,
      uniprot_id  = uniprot_id,
      gene_symbol = gene_symbol
    ),
  chronos_alamar_protein_tibble |>
    dplyr::select(
      platform_id = platform_id,
      uniprot_id  = uniprot_id,
      gene_symbol = gene_symbol
    )
)
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
  for (tbl_name in names(single_protein_tibble_list)) {
    cat("Calculating for single protein table : ", tbl_name, "; dims (", dim(single_protein_tibble_list[[tbl_name]]), ")\n")

    results <- analyze_multiple_prot_marker_by_tibble(
      single_prot_tibble = single_protein_tibble_list[[tbl_name]],
      single_prot_metadata = metadata,
      protein_intensity_log2 = protein_intensity_log2,
      annotations = annotations,
      calc_adj_p = FALSE, ## Adj p vals not used for further analysis -- set to FALSE
      sort = TRUE,
      mapping_cross_olink_platforms = TRUE,
      mapping_cross_somalogic_platforms = TRUE
    ) |>
      dplyr::mutate(
        protein_pair_tbl_name = !!tbl_name
      ) |>
      dplyr::relocate(protein_pair_tbl_name)

    # Collect
    all_results[[tbl_name]] <- results
    cat("Done for single protein table : ", tbl_name, " ; results dims (", dim(results), ")\n")
  }

  output_file_path <- glue::glue("{output_dir}/cohort_{selected_contributor_code}{visit_str}_outliers_rm_{matching_ratio_str}_matched_grp_3tbl_single_protein_validation_results.csv")
  
  final_results <- do.call(rbind, all_results)

  # filter out unmatched or unconverged rows (NA rows)
  final_results <- final_results |>
    dplyr::filter(!is.na(r_squared))

  final_results |>
    readr::write_csv(output_file_path)
  cat("Done for GNPC cohort : ", selected_contributor_code, " ; output dims (", dim(final_results), ")\n")
}
