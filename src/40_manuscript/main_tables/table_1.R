## SPDX-FileCopyrightText: 2026 Alkahest, Inc. <chronospd_data_access@alkahest.com>
##
## SPDX-License-Identifier: MIT

here::i_am("src/40_manuscript/main_tables/table_1.R")

infile <- glue::glue("{get_config()$input_dir}/name=Chronos_phase1_critical_variables_table/version=29432/dataset.parquet")
outfile_donor_characteristics <- glue::glue("{get_config()$output_dir}/main_tables/table_1/donor-level-characteristics.docx")
outfile_sample_characteristics <- glue::glue("{get_config()$output_dir}/main_tables/table_1/sample-level-characteristics.docx")

dir.create(dirname(outfile_donor_characteristics), showWarnings = FALSE, recursive = TRUE)
dir.create(dirname(outfile_sample_characteristics), showWarnings = FALSE, recursive = TRUE)

critical_variables <- infile |>
  arrow::open_dataset() |>
  dplyr::select(
    gender,
    race,
    hypertension,
    diabetes,
    age_at_baseline,
    age_pd_onset,
    group,
    patient_id,
    age_at_donation,
    years_since_pd_onset,
    recent_plasma_donations,
    cluster,
    group,
    sample_id
  ) |>
  dplyr::collect() |>
  dplyr::mutate(
    group = factor(
      group,
      levels = c("pd", "control"),
      labels = c("PD cases", "Controls")
    )
  )

donor_level_characteristics <- critical_variables |>
  dplyr::select(
    gender,
    race,
    hypertension,
    diabetes,
    age_at_baseline,
    age_pd_onset,
    group,
    patient_id
  ) |>
  dplyr::distinct() |>
  dplyr::mutate(
    age_pd_onset = ifelse(group == "PD cases", age_pd_onset, NA)
  ) |>
  gtsummary::tbl_summary(
    include = c(
      gender,
      race,
      hypertension,
      diabetes,
      age_at_baseline,
      age_pd_onset,
      group
    ),
    type = list(c(age_at_baseline, age_pd_onset) ~ "continuous"),
    digits = list(c(gtsummary::all_continuous()) ~ 1),
    statistic = list(c(gtsummary::all_continuous()) ~ "{mean} ({sd},{min},{max}) {median} ({p25},{p75})"),
    by = c("group")
  ) |>
  gtsummary::modify_header(label = "**Donor-level characteristics**")

donor_level_characteristics |>
  gtsummary::as_gt() |>
  gt::gtsave(outfile_donor_characteristics)


samples_per_patient <- critical_variables |>
  dplyr::select(patient_id, sample_id, group) |>
  dplyr::distinct() |>
  dplyr::group_by(patient_id, group) |>
  dplyr::summarise(
    "samples per individual" = dplyr::n(),
    .groups = "drop"
  )
tbl_samples <- samples_per_patient |>
  gtsummary::tbl_summary(
    include = "samples per individual",
    by = group,
    type = "samples per individual" ~ "continuous",
    digits = "samples per individual" ~ 1,
    statistic = "samples per individual" ~ "{mean} ({sd},{min},{max}) {median} ({p25},{p75})"
  )


sample_level_characteristics <- critical_variables |>
  dplyr::select(
    age_at_donation,
    years_since_pd_onset,
    recent_plasma_donations,
    PC_cluster = cluster,
    group,
    sample_id,
    patient_id
  ) |>
  dplyr::mutate(
    years_since_pd_onset = ifelse(group == "PD cases", years_since_pd_onset, NA),
    PC_cluster = dplyr::case_when(
      PC_cluster == 1 ~ "cluster1",
      PC_cluster == 2 ~ "cluster2",
      .default = PC_cluster
    )
  ) |>
  gtsummary::tbl_summary(
    include = c(age_at_donation, years_since_pd_onset, recent_plasma_donations, PC_cluster, group),
    type = list(c(age_at_donation, years_since_pd_onset) ~ "continuous"),
    digits = list(c(gtsummary::all_continuous()) ~ 1),
    statistic = list(c(gtsummary::all_continuous()) ~ "{mean} ({sd},{min},{max}) {median} ({p25},{p75})"),
    by = c("group")
  ) |>
  gtsummary::modify_header(label = "**Sample-level characteristics**")

final_table <- gtsummary::tbl_stack(
  list(sample_level_characteristics, tbl_samples),
  quiet = TRUE
)

final_table |>
  gtsummary::as_gt() |>
  gt::gtsave(outfile_sample_characteristics)
