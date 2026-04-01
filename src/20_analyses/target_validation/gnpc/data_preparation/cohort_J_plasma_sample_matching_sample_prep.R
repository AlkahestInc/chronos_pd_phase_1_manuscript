## SPDX-FileCopyrightText: 2026 Alkahest, Inc. <chronospd_data_access@alkahest.com>
##
## SPDX-License-Identifier: MIT

here::i_am("src/20_analyses/target_validation/gnpc/data_preparation/cohort_J_plasma_sample_matching_sample_prep.R")

input_dir <- "/home/rstudio"

# check or create output path
output_dir <- glue::glue("{get_config()$checkpoint_dir}/gnpc_cohort_cleanup_matched_data")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

annotations <- glue::glue(
  "{input_dir}/data_prep/annotation/annotation_tbl_v6.rds",
) |>
  readr::read_rds()

metadata <- glue::glue(
  "{input_dir}/data_prep/metadata/clinical_metadata_tbl_v6_contributor_J.rds",
) |>
  readr::read_rds() |>
  dplyr::filter(
    sample_matrix == "EDTA Plasma",
    !is.na(age_at_visit),
    !is.na(sex),
    !is.na(ad)
  ) |>
  # Get the sample from the last visit
  dplyr::group_by(person_id) |>
  dplyr::slice_max(
    order_by = sequential_visit_number,
    n = 1,
    with_ties = FALSE
  ) |>
    dplyr::ungroup()

protein_intensity_log2 <- glue::glue(
  "{input_dir}/data_prep/proteomics/proteomic_tbl_v6_contributor_J.rds",
) |>
  readr::read_rds() |>
  dplyr::semi_join(
    metadata,
    by = "sample_id"
  ) |>
  dplyr::mutate(
    dplyr::across(where(is.numeric), ~ log2(.x + 1))
  )

#####
# Detect samples that are PCA outliers
#####
pca_outlier_row_numbers <- mt::pca.outlier(
  protein_intensity_log2 |>
    dplyr::select(-sample_id)
) |>
  dplyr::pull(outlier)

pca_outlier_sample_ids <- protein_intensity_log2 |>
  dplyr::slice(pca_outlier_row_numbers)
cat("Number of PCA outliers :", nrow(pca_outlier_sample_ids), "\n")

#####
# Remove PCA outliers
#####
filtered_protein_intensity_log2 <- protein_intensity_log2 |>
  dplyr::anti_join(
    pca_outlier_sample_ids,
    by = "sample_id"
  )
filtered_metadata <- metadata |>
  dplyr::anti_join(
    pca_outlier_sample_ids,
    by = "sample_id"
  )

metadata_output_stem <- glue::glue("{get_config()$checkpoint_dir}/data_prep/gnpc_cohort_cleanup_matched_data/cohort_J_plasma_last_visit_outliers_rm_1to1_matched_grp_filtered_metadata")
protein_intensity_log2_output_stem <- glue::glue("{get_config()$checkpoint_dir}/data_prep/gnpc_cohort_cleanup_matched_data/cohort_J_plasma_last_visit_outliers_rm_1to1_matched_grp_filtered_protein_intensity_log2")

matched_output <- match_stratified_cohort(
  data = filtered_metadata,
  predictors = c("age_at_visit", "sex", "ad"), ## note here the predictors and response is for sample matching
  response = "pd",
  tolerance_in_age = 5,
  ratio = 1,
  seed = 42,
  method = "nearest",
  distance = "mahalanobis"
)

# Run regression on matched dataset
validate_lm <- stats::glm(
  pd ~ age_at_visit + sex,
  data = matched_output$matched_metadata,
  family = binomial()
)

matched_filtered_metadata <- matched_output$matched_metadata
cat("Number of matched samples :", nrow(matched_filtered_metadata), "\n")

matched_filtered_protein_intensity_log2 <- filtered_protein_intensity_log2 |>
dplyr::semi_join(
  matched_filtered_metadata,
  by = "sample_id"
)

matched_filtered_protein_intensity_log2 |>
readr::write_csv(glue::glue("{protein_intensity_log2_output_stem}.csv"))
matched_filtered_protein_intensity_log2 |>
readr::write_rds(glue::glue("{protein_intensity_log2_output_stem}.rds"))

matched_filtered_metadata |>
readr::write_csv(glue::glue("{metadata_output_stem}.csv"))
matched_filtered_metadata |>
readr::write_rds(glue::glue("{metadata_output_stem}.rds"))
