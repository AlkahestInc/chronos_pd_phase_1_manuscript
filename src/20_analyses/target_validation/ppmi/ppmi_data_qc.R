## SPDX-FileCopyrightText: 2026 Alkahest, Inc. <chronospd_data_access@alkahest.com>
##
## SPDX-License-Identifier: MIT

here::i_am("src/20_analyses/target_validation/ppmi/ppmi_data_qc.R")

# Load the preprocessed PPMI data
p9k <- glue::glue(
  "{get_config()$checkpoint_dir}/name={name}/version={version}/dataset.rds",
  name = "PPMI_p9000_Plasma_MS/",
  version = "29985"
) |>
  readr::read_rds()

# patch: add ID column to dataProt
p9k$dataProt$ID <- rownames(p9k$dataProt)

# patch: add protein column to nom
p9k$nom <- p9k$nom |>
  dplyr::mutate(protein = paste(LABEL, "NPX", sep = "_"))

# patch: add age_at_first_visit variable to donor table
p9k$pheno.don$age_at_first_visit <- NA
for (k in seq_len(nrow(p9k$pheno.don))) {
  i <- which(p9k$pheno.samp$PATNO == p9k$pheno.don$PATNO[k])
  p9k$pheno.don$age_at_first_visit[k] <- min(p9k$pheno.samp$age_at_visit[i], na.rm = TRUE)
}

# generate LD dataset, all visits ####

# baseline meta using earliest date per subject ####
PPMI_pd_meta_BL <- p9k$pheno.samp |>
  dplyr::group_by(PATNO) |>
  dplyr::filter(
    cohort %in% c("control", "PD"),
    visit_date == min(visit_date)
  ) |> 
  dplyr::left_join(p9k$pheno.don |>
    dplyr::select(
      PATNO, ENROLL_STATUS, sex, race, education_years, smoking, alcohol, handedness,
      hypertension, diabetes, lipid_cholesterol, depression_anxiety, cognitive_deficit,
      pddiag, age_at_pd_diag, early_onset_pd, age_at_first_visit
    ), by = "PATNO") |>
  dplyr::rename(ID = SampleId) |>
  dplyr::ungroup()

# LD meta
PPMI_pd_meta_LD <- p9k$pheno.samp |>
  dplyr::filter(cohort %in% c("control", "PD")) |>
  dplyr::group_by(PATNO) |>
  dplyr::mutate(bl_subjects = ifelse(any(EVENT_ID == "BL"), "with_bl", "no_bl")) |>
  dplyr::ungroup() |>
  dplyr::left_join(p9k$pheno.don |>
    dplyr::select(
      PATNO, ENROLL_STATUS, sex, race, education_years, smoking, alcohol, handedness,
      hypertension, diabetes, lipid_cholesterol, depression_anxiety, cognitive_deficit,
      pddiag, age_at_pd_diag, early_onset_pd, age_at_first_visit
    ), by = "PATNO") |>
  dplyr::rename(ID = SampleId) |>
  dplyr::ungroup()

# LD proteomics
PPMI_pd_protein_LD <- p9k$dataProt |>
  dplyr::filter(ID %in% PPMI_pd_meta_LD$ID)
rownames(PPMI_pd_protein_LD) <- NULL

# BL proteomics
PPMI_pd_protein_BL <- PPMI_pd_protein_LD |>
  dplyr::filter(ID %in% PPMI_pd_meta_BL$ID)

# merged LD dataset ####
merged_data_LD <- PPMI_pd_protein_LD |>
  dplyr::inner_join(PPMI_pd_meta_LD, by = "ID")

# protein names
protein_names <- paste0(p9k$nom$LABEL, "_NPX")

# QC steps ####

# missingness ####
# Identify variable types
numeric_vars_LD <- PPMI_pd_meta_LD |>
  dplyr::select(where(is.numeric)) |>
  names()

categorical_vars_LD <- PPMI_pd_meta_LD |>
  dplyr::select(where(~ is.character(.) | is.factor(.))) |>
  names()

# Missingness summary for numeric variables
missing_numeric_LD <- PPMI_pd_meta_LD |>
  dplyr::summarise(across(all_of(numeric_vars_LD), ~ mean(is.na(.)))) |>
  tidyr::pivot_longer(
    cols = everything(),
    names_to = "variable",
    values_to = "missing_pct"
  ) |>
  dplyr::mutate(type = "numeric")

# Missingness summary for categorical variables
missing_categorical_LD <- PPMI_pd_meta_LD |>
  dplyr::summarise(across(all_of(categorical_vars_LD), ~ mean(is.na(.)))) |>
  tidyr::pivot_longer(
    cols = everything(),
    names_to = "variable",
    values_to = "missing_pct"
  ) |>
  dplyr::mutate(type = "categorical")

# Combine both summaries
missing_summary_LD <- dplyr::bind_rows(missing_numeric_LD, missing_categorical_LD)

# View
ggplot2::ggplot(
  missing_summary_LD |> dplyr::filter(missing_pct > 0.0),
  ggplot2::aes(x = reorder(variable, missing_pct), y = missing_pct, fill = type)
) +
  ggplot2::geom_col() +
  ggplot2::coord_flip() +
  ggplot2::scale_y_continuous(labels = scales::percent) +
  ggplot2::labs(
    title = "Column Missingness by Variable Type",
    x = "Variable",
    y = "Missingness (%)"
  ) +
  ggplot2::theme_minimal()

# do not need to summarize variables with >10% missingness
# subject level table 1 ####
summary_cols_CS_BL <- c(
  "cohort", "age_at_visit", "sex", "race", "education_years", "BMI", "smoking", "alcohol",
  "handedness", "hypertension", "diabetes", "lipid_cholesterol", "depression_anxiety", "cognitive_deficit",
  "age_at_pd_diag", "early_onset_pd", "age_at_first_visit"
)
continuous_cols_CS_BL <- c("age_at_visit", "education_years", "BMI", "age_at_pd_diag", "age_at_first_visit")

PPMI_pd_meta_BL |>
  dplyr::select(all_of(summary_cols_CS_BL)) |>
  dplyr::distinct() |>
  gtsummary::tbl_summary(
    include = c(all_of(summary_cols_CS_BL)),
    type = list(continuous_cols_CS_BL ~ "continuous"),
    digits = list(c(gtsummary::all_continuous()) ~ 1),
    statistic = list(c(gtsummary::all_continuous()) ~ "{mean} ({sd},{min},{max}) {median} ({p25},{p75})"),
    by = c("cohort")
  )


# subject level table 1 ####
summary_cols_LD <- c(
  "cohort", "EVENT_ID", "age_at_visit", "sex", "race", "education_years", "BMI", "smoking", "alcohol",
  "handedness", "hypertension", "diabetes", "lipid_cholesterol", "depression_anxiety", "cognitive_deficit",
  "age_at_pd_diag", "early_onset_pd", "age_at_first_visit"
)
continuous_cols_LD <- c("age_at_visit", "education_years", "BMI", "age_at_pd_diag", "age_at_first_visit")

PPMI_pd_meta_LD |>
  dplyr::select(
    tidyselect::all_of(summary_cols_LD)
  ) |>
  dplyr::distinct() |>
  gtsummary::tbl_summary(
    include = tidyselect::all_of(summary_cols_LD),
    type = list(continuous_cols_LD ~ "continuous"),
    digits = list(c(gtsummary::all_continuous()) ~ 1),
    statistic = list(c(gtsummary::all_continuous()) ~ "{mean} ({sd},{min},{max}) {median} ({p25},{p75})"),
    by = c("cohort")
  )

# PCA for outlier and local structures in proteomics ####
# Calculate number of missing values per protein ####
proteomic_missing <- PPMI_pd_protein_LD |>
  dplyr::select(-ID) |> # exclude ID column if present
  dplyr::summarise(across(everything(), ~ sum(is.na(.)))) |>
  tidyr::pivot_longer(cols = everything(), names_to = "variable", values_to = "missing_count")

# remove proteins with missingness
PPMI_pd_protein_LD[1:3, 1460:1464]
remove_index <- vector()
for (i in 1:1463) {
  if (sum(is.na(PPMI_pd_protein_LD[, i])) >= 1) {
    remove_index <- c(remove_index, i)
    print(paste(i, sum(is.na(PPMI_pd_protein_LD[, i]))))
  }
}
remove_index
for (i in remove_index) {
  print(paste(i, which(is.na(PPMI_pd_protein_LD[, i]))))
}
not_protein_name_index <- c(remove_index, 1464) # 1464 for ID index
protein_names_NA_removed <- colnames(PPMI_pd_protein_LD[, -not_protein_name_index]) # here I removed all proteins with at least one missing sample, maybe just remove the last two with 139 missing samples
length(protein_names_NA_removed) # 1456

Pro_PCA_meta1(
  dataset = merged_data_LD,
  protein_names = protein_names_NA_removed,
  var_names = "cohort",
  title_name = paste("PPMI", "cohort", sep = ""),
  cont = "cate",
  nsd = 6,
  IDs_for_outliers = c("ID", "PATNO", "cohort", "age_at_visit")
)

# based on sample median for outliers ####
sample_median <- apply(merged_data_LD |>
  dplyr::select(all_of(protein_names_NA_removed)), 1, median)
mean_protein <- mean(sample_median)
sd_protein <- sd(sample_median)

plot_data <- data.frame(
  sample_median = sample_median,
  group = merged_data_LD$cohort,
  EVENT_ID = merged_data_LD$EVENT_ID,
  outlier_sd = abs(sample_median - mean_protein) > 5 * sd_protein
)
# check outliers
# no outlier detected
# ploat sample median per plate
ggplot2::ggplot(plot_data, ggplot2::aes(x = EVENT_ID, y = sample_median)) +
  ggplot2::geom_jitter(width = 0.3, ggplot2::aes(color = outlier_sd)) +
  ggplot2::geom_boxplot(ggplot2::aes(fill = factor(EVENT_ID)), outlier.shape = NA) +
  ggplot2::geom_hline(yintercept = mean_protein, lty = 2, col = "blue") +
  ggplot2::theme_classic()


# summary
table(merged_data_LD$cohort, merged_data_LD$EVENT_ID)
length(unique(merged_data_LD$PATNO[which(merged_data_LD$cohort == "PD")]))
analysis_ready_table <- merged_data_LD |> dplyr::filter(EVENT_ID %in% c("BL", "V02", "V04", "V06", "V10")) # number of samples>=20
table(analysis_ready_table$cohort)
length(unique(analysis_ready_table$PATNO[which(analysis_ready_table$cohort == "PD")]))
length(unique(analysis_ready_table$PATNO[which(analysis_ready_table$cohort == "control")]))
