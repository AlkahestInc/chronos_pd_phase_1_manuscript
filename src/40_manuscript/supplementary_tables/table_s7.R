## SPDX-FileCopyrightText: 2026 Alkahest, Inc. <chronospd_data_access@alkahest.com>
##
## SPDX-License-Identifier: MIT

here::i_am("src/40_manuscript/supplementary_tables/table_s7.R")

library(lspline)

# ==============================================================================
# Input file paths
# ==============================================================================
# Analysis-ready tibbles
somalogic_tibble_path <- glue::glue("{get_config()$input_dir}/name=Chronos_phase1_proteomics_somalogic_analysis_ready_three_tibble/version=29650/dataset.rds")
olink_tibble_path <- glue::glue("{get_config()$input_dir}/name=Chronos_phase1_proteomics_olink_analysis_ready_three_tibble/version=29465/dataset.rds")
alamar_tibble_path <- glue::glue("{get_config()$input_dir}/name=Chronos_phase1_proteomics_alamar_analysis_ready_three_tibble/version=29649/dataset.rds")

# SomaScan data paths
soma_all_path <- glue::glue("{get_config()$input_dir}/name=chronos_somalogic_all-pd_master_summary/version=29806/dataset.parquet")
soma_pre_path <- glue::glue("{get_config()$input_dir}/name=chronos_somalogic_pre_onset-pd_master_summary/version=29807/dataset.parquet")
soma_post_path <- glue::glue("{get_config()$input_dir}/name=chronos_somalogic_post_onset-pd_master_summary/version=29808/dataset.parquet")
soma_jm_path <- glue::glue("{get_config()$input_dir}/name=chronos_somalogic_joint_model-pd_master_summary/version=29819/dataset.parquet")
soma_model_path <- glue::glue("{get_config()$input_dir}/name=somalogic_screen_model_comparison/version=none/dataset.rds")


# Olink data paths
olink_all_path <- glue::glue("{get_config()$input_dir}/name=chronos_olink_all-pd_master_summary/version=29809/dataset.parquet")
olink_pre_path <- glue::glue("{get_config()$input_dir}/name=chronos_olink_all_pre_onset-pd_master_summary/version=29810/dataset.parquet")
olink_post_path <- glue::glue("{get_config()$input_dir}/name=chronos_olink_all_post_onset-pd_master_summary/version=29811/dataset.parquet")
olink_jm_path <- glue::glue("{get_config()$input_dir}/name=chronos_olink_joint_model-pd_master_summary/version=29820/dataset.parquet")
olink_model_path <- glue::glue("{get_config()$input_dir}/name=olink_screen_model_comparison/version=none/dataset.rds")

# AlamarBio data paths
alamar_all_path <- glue::glue("{get_config()$input_dir}/name=chronos_alamar_all-pd_master_summary/version=29812/dataset.parquet")
alamar_pre_path <- glue::glue("{get_config()$input_dir}/name=chronos_alamar_pre_onset-pd_master_summary/version=29813/dataset.parquet")
alamar_post_path <- glue::glue("{get_config()$input_dir}/name=chronos_alamar_post_onset-pd_master_summary/version=29814/dataset.parquet")
alamar_jm_path <- glue::glue("{get_config()$input_dir}/name=chronos_alamar_joint_model-pd_master_summary/version=29821/dataset.parquet")
alamar_model_path <- glue::glue("{get_config()$input_dir}/name=alamar_screen_model_comparison/version=none/dataset.rds")

# Nonlinear mixed model results paths
somalogic_all_ls_path <- glue::glue("{get_config()$input_dir}/name=somalogic_all_model_results/version=none/dataset.rds")
olink_all_ls_path <- glue::glue("{get_config()$input_dir}/name=olink_all_model_results/version=none/dataset.rds")
alamar_all_ls_path <- glue::glue("{get_config()$input_dir}/name=alamar_all_model_results/version=none/dataset.rds")

Olink_med_table_path <- glue::glue("{get_config()$input_dir}/name=Olink_Levodopa_comparison_Phase1_PD_MSvers_LMMs_simplified/version=29945/dataset.parquet")
SomaLogic_med_table_path <- glue::glue("{get_config()$input_dir}/name=Soma_Levodopa_comparison_Phase1_PD_MSvers_LMMs_simplified/version=29944/dataset.parquet")
Alamar_med_table_path <- glue::glue("{get_config()$input_dir}/name=AlamarBio_Levodopa_comparison_Phase1_PD_MSvers_LMMs_simplified/version=29946/dataset.parquet")

# ==============================================================================
# Load input data
# ==============================================================================
# Analysis-ready tibbles
somalogic_tibble <- readr::read_rds(somalogic_tibble_path)
olink_tibble <- readr::read_rds(olink_tibble_path)
alamar_tibble <- readr::read_rds(alamar_tibble_path)

# SomaScan data
soma_all <- arrow::read_parquet(soma_all_path)
soma_pre <- arrow::read_parquet(soma_pre_path)
soma_post <- arrow::read_parquet(soma_post_path)
soma_jm <- arrow::read_parquet(soma_jm_path)
soma_model <- readr::read_rds(soma_model_path)

# Olink data
olink_all <- arrow::read_parquet(olink_all_path)
olink_pre <- arrow::read_parquet(olink_pre_path)
olink_post <- arrow::read_parquet(olink_post_path)
olink_jm <- arrow::read_parquet(olink_jm_path)
olink_model <- readr::read_rds(olink_model_path)

# AlamarBio data
alamar_all <- arrow::read_parquet(alamar_all_path)
alamar_pre <- arrow::read_parquet(alamar_pre_path)
alamar_post <- arrow::read_parquet(alamar_post_path)
alamar_jm <- arrow::read_parquet(alamar_jm_path)
alamar_model <- readr::read_rds(alamar_model_path)

# Nonlinear mixed model results
somalogic_all_ls <- readr::read_rds(somalogic_all_ls_path)
olink_all_ls <- readr::read_rds(olink_all_ls_path)
alamar_all_ls <- readr::read_rds(alamar_all_ls_path)

# ==============================================================================
# Define output files and helper functions
# ==============================================================================
output_file <- glue::glue("{get_config()$output_dir}/supplementary_tables/table_s7.tsv")
dir.create(dirname(output_file), showWarnings = FALSE, recursive = TRUE)

"%nin%" <- Negate("%in%")

# ==============================================================================
# Prepare data tables
# ==============================================================================

somalogic_table <- somalogic_tibble$meta_data |>
  dplyr::inner_join(
    somalogic_tibble$prot_data,
    by = "ID"
  )

olink_table <- olink_tibble$meta_data |>
  dplyr::inner_join(
    olink_tibble$prot_data,
    by = "ID"
  )

alamar_table <- alamar_tibble$meta_data |>
  dplyr::inner_join(
    alamar_tibble$prot_data,
    by = "ID"
  )

# ==============================================================================
# Select proteins with significant BIC differences for nonlinear models
# ==============================================================================

somalogic_selected_proteins_BIC <- names(somalogic_all_ls)[
  sapply(somalogic_all_ls, function(df) any(df$diff_lsBIC > 2, na.rm = TRUE))
] # 53

somalogic_all_ls_selected_BIC <- dplyr::bind_rows(somalogic_all_ls[somalogic_selected_proteins_BIC], .id = "protein") |>
  dplyr::filter(diff_lsBIC > 2 & group_diff_sig != "neither") |>
  dplyr::group_by(protein) |>
  dplyr::filter(diff_lsBIC == max(diff_lsBIC))

# Olink
olink_selected_proteins_BIC <- names(olink_all_ls)[
  sapply(olink_all_ls, function(df) any(df$diff_lsBIC > 2, na.rm = TRUE))
]
olink_all_ls_selected_BIC <- dplyr::bind_rows(olink_all_ls[olink_selected_proteins_BIC], .id = "protein") |>
  dplyr::filter(diff_lsBIC > 2 & group_diff_sig != "neither") |>
  dplyr::group_by(protein) |>
  dplyr::filter(diff_lsBIC == max(diff_lsBIC))

olink_all_ls_selected_BIC # one

### alamarbio ####
alamar_selected_proteins_BIC <- names(alamar_all_ls)[
  sapply(alamar_all_ls, function(df) any(df$diff_lsBIC > 2, na.rm = TRUE))
]

alamar_all_ls_selected_BIC <- dplyr::bind_rows(alamar_all_ls[alamar_selected_proteins_BIC], .id = "protein") |>
  dplyr::filter(diff_lsBIC > 2 & group_diff_sig != "neither") |>
  dplyr::group_by(protein) |>
  dplyr::filter(diff_lsBIC == max(diff_lsBIC)) # 2 with significant slope difference, ARSA (-1),  DDC(3)
alamar_all_ls_selected_BIC

# ==============================================================================
# Create cross-platform combined results
# ==============================================================================
all_combined <- soma_all |>
  dplyr::select(probe_id, estimate, cohens_d, p_val, q_val, uniprot_id, gene_symbol) |>
  dplyr::distinct() |>
  dplyr::mutate(platform = "SomaScanV5") |>
  dplyr::bind_rows(
    olink_all |>
      dplyr::select(probe_id, estimate, cohens_d, p_val, q_val, uniprot_id, gene_symbol) |>
      dplyr::distinct() |>
      dplyr::mutate(platform = "Olink_HT")
  ) |>
  dplyr::bind_rows(
    alamar_all |>
      dplyr::select(probe_id, estimate, cohens_d, p_val, q_val, uniprot_id, gene_symbol) |>
      dplyr::distinct() |>
      dplyr::mutate(platform = "Alamarbio_CNS")
  )


### mixed model from preDx samples
preDx_combined <- soma_pre |>
  dplyr::select(probe_id, estimate, cohens_d, p_val, q_val, uniprot_id, gene_symbol) |>
  dplyr::distinct() |>
  dplyr::mutate(platform = "SomaScanV5") |>
  dplyr::bind_rows(olink_pre |>
    dplyr::select(probe_id, estimate, cohens_d, p_val, q_val, uniprot_id, gene_symbol) |>
    dplyr::distinct() |>
    dplyr::mutate(platform = "Olink_HT")) |>
  dplyr::bind_rows(alamar_pre |>
    dplyr::select(probe_id, estimate, cohens_d, p_val, q_val, uniprot_id, gene_symbol) |>
    dplyr::distinct() |>
    dplyr::mutate(platform = "Alamarbio_CNS"))

### mixed model from postDx samples
postDx_combined <- soma_post |>
  dplyr::select(probe_id, estimate, cohens_d, p_val, q_val, uniprot_id, gene_symbol) |>
  dplyr::distinct() |>
  dplyr::mutate(platform = "SomaScanV5") |>
  dplyr::bind_rows(olink_post |>
    dplyr::select(probe_id, estimate, cohens_d, p_val, q_val, uniprot_id, gene_symbol) |>
    dplyr::distinct() |>
    dplyr::mutate(platform = "Olink_HT")) |>
  dplyr::bind_rows(alamar_post |>
    dplyr::select(probe_id, estimate, cohens_d, p_val, q_val, uniprot_id, gene_symbol) |>
    dplyr::distinct() |>
    dplyr::mutate(platform = "Alamarbio_CNS"))

### joint model results
risk_combined <- soma_jm |>
  dplyr::filter(rhat < 1.05) |>
  dplyr::select(probe_id, uniprot_id, gene_symbol, estimate, p_val, q_val, risk, risk_lower, risk_upper) |>
  dplyr::mutate(platform = "SomaScanV5") |>
  dplyr::bind_rows(olink_jm |>
    dplyr::filter(rhat < 1.05) |>
    dplyr::select(probe_id, uniprot_id, gene_symbol, estimate, p_val, q_val, risk, risk_lower, risk_upper) |>
    dplyr::mutate(platform = "Olink_HT")) |>
  dplyr::bind_rows(alamar_jm |>
    dplyr::filter(rhat < 1.05) |>
    dplyr::select(probe_id, uniprot_id, gene_symbol, estimate, p_val, q_val, risk, risk_lower, risk_upper) |>
    dplyr::mutate(platform = "Alamarbio_CNS"))


### linear spline model results combined
linear_spline_model_selection <- somalogic_all_ls_selected_BIC |>
  dplyr::mutate(platform = "SomaScanV5") |>
  dplyr::bind_rows(olink_all_ls_selected_BIC |>
    dplyr::mutate(platform = "Olink_HT")) |>
  dplyr::bind_rows(alamar_all_ls_selected_BIC |>
    dplyr::mutate(platform = "Alamarbio_CNS")) |>
  dplyr::rename(probe_id = protein) |>
  dplyr::inner_join(all_combined |>
    dplyr::select(probe_id, uniprot_id, gene_symbol), by = "probe_id")


# ==============================================================================
# Discovery of PD biomarker trajectories
# ==============================================================================
# 5 groups identified and then reduced to 3 groups as reported in manuscripts
# Based on linear mixed model results and non-linear model results

# Group 1: Transient preclinical biomarkers
### Transient preclinical biomarkers:
## Identified PD biomarkers with significantly different rate of change in the first slope
## with the inflection point before PD onset or around PD onset from linear spline mixed model results

transient_preclinical_soma <- somalogic_all_ls_selected_BIC |>
  dplyr::filter(Knot <= 0 & group_diff_sig == "first")
transient_preclinical_olink <- olink_all_ls_selected_BIC |>
  dplyr::filter(Knot <= 0 & group_diff_sig == "first")
transient_preclinical_alamar <- alamar_all_ls_selected_BIC |>
  dplyr::filter(Knot <= 0 & group_diff_sig == "first")

# List of transient_preclinical_probe_id from all 3 platforms
transient_preclinical_probe_id <- c(
  transient_preclinical_soma$protein,
  transient_preclinical_olink$protein,
  transient_preclinical_alamar$protein
)

# ------------------------------------------------------------------------------
# Refit spline models for transient preclinical biomarkers
# ------------------------------------------------------------------------------
### Somalogic
# Set custom control parameters and do the spline modeling
ctrl <- lme4::lmerControl(
  optimizer = "bobyqa",
  optCtrl = list(maxfun = 100000)
)

spline_results_soma <- data.frame()
predicted_soma <- data.frame()
for (i in 1:nrow(transient_preclinical_soma)) { #
  # modeling data table generation
  my_table <- somalogic_table |>
    dplyr::select(
      all_of(transient_preclinical_soma$protein[i]),
      years_since_pd_onset, group, patient_id, age_at_baseline, gender, recent_plasma_donations, cluster
    )
  names(my_table)[1] <- "protein"
  # fit linear spline model with selected knot
  lmer.linear_spline1 <- lmerTest::lmer(
    protein ~ age_at_baseline + gender + recent_plasma_donations + cluster +
      lspline::lspline(years_since_pd_onset, transient_preclinical_soma$Knot[i]) * group +
      (1 | patient_id),
    REML = T, data = my_table, control = ctrl
  )
  A <- summary(lmer.linear_spline1)
  Ttable <- as.data.frame(A$coefficients)
  Ttable$coef.name <- rownames(Ttable)
  Ttable$protein <- transient_preclinical_soma$protein[i]
  Ttable$inflection_point <- transient_preclinical_soma$Knot[i]
  spline_results_soma <- spline_results_soma |>
    dplyr::bind_rows(Ttable)

  # prediction and normalization
  # Get estimated marginal means and 95% CI
  # Create a prediction grid

  pred_grid <- expand.grid(
    years_since_pd_onset = seq(-12, 9, by = 0.1),
    group = unique(somalogic_table$group),
    age_at_baseline = mean(somalogic_table$age_at_baseline), # Set to reference value or mean value
    gender = "female", # Set to reference value or mean value
    recent_plasma_donations = "at least one donation", # Set to reference value or mean value
    cluster = 1 # Set to reference value or mean value
  ) |>
    dplyr::mutate(cluster = as.factor(cluster))

  # Add spline columns
  pred_grid_splines <- lspline::lspline(pred_grid$years_since_pd_onset,
    knots = transient_preclinical_soma$Knot[i]
  )
  colnames(pred_grid_splines) <- c("spline1", "spline2")

  pred_grid <- cbind(pred_grid, pred_grid_splines)

  pred_grid$pred <- predict(lmer.linear_spline1, newdata = pred_grid, re.form = NA)


  df_wide <- tidyr::pivot_wider(
    pred_grid |>
      dplyr::select(group, years_since_pd_onset, pred) |>
      dplyr::distinct(),
    names_from = group,
    values_from = pred
  )


  df_wide$normalized <- df_wide$pd - df_wide$control

  df_wide$probe_id <- transient_preclinical_soma$protein[i]
  predicted_soma <- predicted_soma |>
    dplyr::bind_rows(df_wide |>
      dplyr::select(years_since_pd_onset, normalized, probe_id))
}
names(spline_results_soma)[c(2, 5, 7)] <- c("SE", "p_value", "probe_id")
names(spline_results_soma)
spline_results_soma_updated <- spline_results_soma |>
  dplyr::filter(coef.name == "lspline::lspline(years_since_pd_onset, transient_preclinical_soma$Knot[i])1:grouppd") |>
  dplyr::select(probe_id, Estimate, p_value, inflection_point) |>
  dplyr::inner_join(soma_all |>
    dplyr::select(probe_id, gene_symbol, uniprot_id), by = "probe_id")
spline_results_soma_updated |>
  dplyr::arrange(probe_id)

predicted_soma <- predicted_soma |>
  dplyr::inner_join(soma_all |>
    dplyr::select(probe_id, gene_symbol, uniprot_id), by = "probe_id")

### Olink
spline_results_olink <- data.frame()
predicted_olink <- data.frame()

for (i in 1:nrow(transient_preclinical_olink)) {
  my_table <- olink_table |>
    dplyr::select(
      all_of(transient_preclinical_olink$protein[i]),
      years_since_pd_onset, group, patient_id, age_at_baseline, gender, recent_plasma_donations, cluster
    )
  names(my_table)[1] <- "protein"
  lmer.linear_spline1 <- lmerTest::lmer(
    protein ~ age_at_baseline + gender + recent_plasma_donations + cluster +
      lspline::lspline(years_since_pd_onset, transient_preclinical_olink$Knot[i]) * group +
      (1 | patient_id),
    REML = T, data = my_table, control = ctrl
  )
  A <- summary(lmer.linear_spline1)
  Ttable <- as.data.frame(A$coefficients)
  Ttable$coef.name <- rownames(Ttable)
  Ttable$protein <- transient_preclinical_olink$protein[i]
  Ttable$inflection_point <- transient_preclinical_olink$Knot[i]
  spline_results_olink <- spline_results_olink |>
    dplyr::bind_rows(Ttable)

  # prediction and normalization
  # Get estimated marginal means and 95% CI
  # Create a prediction grid
  pred_grid <- expand.grid(
    years_since_pd_onset = seq(-12, 9, by = 0.1),
    group = unique(olink_table$group),
    age_at_baseline = mean(olink_table$age_at_baseline), # Set to reference or mean value
    gender = "female",
    recent_plasma_donations = "at least one donation",
    cluster = 1
  ) |>
    dplyr::mutate(cluster = as.factor(cluster))

  # Add spline columns
  pred_grid_splines <- lspline::lspline(pred_grid$years_since_pd_onset,
    knots = transient_preclinical_olink$Knot[i]
  )
  colnames(pred_grid_splines) <- c("spline1", "spline2")

  pred_grid <- cbind(pred_grid, pred_grid_splines)

  pred_grid$pred <- predict(lmer.linear_spline1, newdata = pred_grid, re.form = NA)


  df_wide <- tidyr::pivot_wider(
    pred_grid |>
      dplyr::select(group, years_since_pd_onset, pred) |>
      dplyr::distinct(),
    names_from = group,
    values_from = pred
  )

  df_wide$normalized <- df_wide$pd - df_wide$control

  df_wide$probe_id <- transient_preclinical_olink$protein[i]
  predicted_olink <- predicted_olink |>
    dplyr::bind_rows(df_wide |>
      dplyr::select(years_since_pd_onset, normalized, probe_id))
}
names(spline_results_olink)[c(2, 5, 7)] <- c("SE", "p_value", "probe_id")

spline_results_olink_updated <- spline_results_olink |>
  dplyr::filter(coef.name == "lspline::lspline(years_since_pd_onset, transient_preclinical_olink$Knot[i])1:grouppd") |>
  dplyr::select(probe_id, Estimate, p_value, inflection_point) |>
  dplyr::inner_join(olink_all |>
    dplyr::select(probe_id, gene_symbol, uniprot_id), by = "probe_id")
predicted_olink <- predicted_olink |>
  dplyr::inner_join(olink_all |>
    dplyr::select(probe_id, gene_symbol, uniprot_id), by = "probe_id")

### Alamar
spline_results_alamar <- data.frame()
predicted_alamar <- data.frame()
for (i in 1:nrow(transient_preclinical_alamar)) {
  my_table <- alamar_table |>
    dplyr::select(
      all_of(transient_preclinical_alamar$protein[i]),
      years_since_pd_onset, group, patient_id, age_at_baseline, gender, recent_plasma_donations, cluster
    )
  names(my_table)[1] <- "protein"
  lmer.linear_spline1 <- lmerTest::lmer(
    protein ~ age_at_baseline + gender + recent_plasma_donations + cluster +
      lspline::lspline(years_since_pd_onset, transient_preclinical_alamar$Knot[i]) * group +
      (1 | patient_id),
    REML = T, data = my_table, control = ctrl
  )
  A <- summary(lmer.linear_spline1)
  Ttable <- as.data.frame(A$coefficients)
  Ttable$coef.name <- rownames(Ttable)
  Ttable$protein <- transient_preclinical_alamar$protein[i]
  Ttable$inflection_point <- transient_preclinical_alamar$Knot[i]
  spline_results_alamar <- spline_results_alamar |>
    dplyr::bind_rows(Ttable)

  # prediction and normalization
  # Get estimated marginal means and 95% CI
  # Create a prediction grid
  pred_grid <- expand.grid(
    years_since_pd_onset = seq(-12, 9, by = 0.1),
    group = unique(alamar_table$group),
    age_at_baseline = mean(alamar_table$age_at_baseline), # Set to reference or mean value
    gender = "female",
    recent_plasma_donations = "at least one donation",
    cluster = 1
  ) |>
    dplyr::mutate(cluster = as.factor(cluster))
  str(pred_grid)
  # Add spline columns
  pred_grid_splines <- lspline::lspline(pred_grid$years_since_pd_onset,
    knots = transient_preclinical_alamar$Knot[i]
  )
  colnames(pred_grid_splines) <- c("spline1", "spline2")

  pred_grid <- cbind(pred_grid, pred_grid_splines)
  head(pred_grid)
  pred_grid$pred <- predict(lmer.linear_spline1, newdata = pred_grid, re.form = NA)


  df_wide <- tidyr::pivot_wider(
    pred_grid |>
      dplyr::select(group, years_since_pd_onset, pred) |>
      dplyr::distinct(),
    names_from = group,
    values_from = pred
  )

  df_wide$normalized <- df_wide$pd - df_wide$control

  df_wide$probe_id <- transient_preclinical_alamar$protein[i]
  predicted_alamar <- predicted_alamar |>
    dplyr::bind_rows(df_wide |>
      dplyr::select(years_since_pd_onset, normalized, probe_id))
}
names(spline_results_alamar)[c(2, 5, 7)] <- c("SE", "p_value", "probe_id")

spline_results_alamar_updated <- spline_results_alamar |>
  dplyr::filter(coef.name == "lspline::lspline(years_since_pd_onset, transient_preclinical_alamar$Knot[i])1:grouppd") |>
  dplyr::select(probe_id, Estimate, p_value, inflection_point) |>
  dplyr::inner_join(alamar_all |>
    dplyr::select(probe_id, gene_symbol, uniprot_id), by = "probe_id")
spline_results_alamar_updated |>
  dplyr::arrange(probe_id)

predicted_alamar <- predicted_alamar |>
  dplyr::inner_join(alamar_all |>
    dplyr::select(probe_id, gene_symbol, uniprot_id), by = "probe_id")

# ------------------------------------------------------------------------------
# Save combined results for transient preclinical biomarkers
# ------------------------------------------------------------------------------
# Mixed model results
all_mixed_model_results <- all_combined |>
  dplyr::filter(probe_id %in% transient_preclinical_probe_id)
all_mixed_model_results$uniprot_id[duplicated(all_mixed_model_results$uniprot_id)] # all unique proteins
# linear spline model results
spline_results_combined <- spline_results_soma_updated |>
  dplyr::bind_rows(spline_results_olink_updated) |>
  dplyr::bind_rows(spline_results_alamar_updated) |>
  dplyr::rename(slope1_group_difference = Estimate, slope1_difference_p_val = p_value)
# combined results
results_transient_preclinical <- all_mixed_model_results |>
  dplyr::select(platform, probe_id, uniprot_id, gene_symbol, estimate, cohens_d, p_val, q_val) |>
  dplyr::rename(
    log2_FC_all = estimate,
    cohens_d_all = cohens_d,
    p_val_all = p_val,
    q_val_all = q_val
  ) |>
  dplyr::inner_join(spline_results_combined |> dplyr::select(-uniprot_id, -gene_symbol), by = "probe_id")


# ==============================================================================
# Group 2: Transitional biomarkers
# ==============================================================================
### Based on linear mixed model results and non-linear models
### Transitional biomarkers:
## Identified PD biomarkers with significantly different rate of change in the second slope
## with inflection point before PD onset or around PD onset from linear spline mixed model results

### Only Somalogic with candidate transitional biomarkers based on the definition
transitional_soma <- somalogic_all_ls_selected_BIC |>
  dplyr::filter(Knot <= 0 & group_diff_sig == "second")
transitional_soma

transitional_probe_id <- transitional_soma$protein

### Set custom control parameters and do the spline modeling ####
ctrl <- lme4::lmerControl(
  optimizer = "bobyqa",
  optCtrl = list(maxfun = 100000)
)

spline_results_soma2 <- data.frame()
predicted_soma2 <- data.frame()
for (i in 1:nrow(transitional_soma)) {
  my_table <- somalogic_table |>
    dplyr::select(
      all_of(transitional_soma$protein[i]),
      years_since_pd_onset, group, patient_id, age_at_baseline, gender, recent_plasma_donations, cluster
    )
  names(my_table)[1] <- "protein"
  lmer.linear_spline1 <- lmerTest::lmer(
    protein ~ age_at_baseline + gender + recent_plasma_donations + cluster +
      lspline::lspline(years_since_pd_onset, transitional_soma$Knot[i]) * group + (1 | patient_id),
    REML = T, data = my_table, control = ctrl
  )
  A <- summary(lmer.linear_spline1)
  Ttable <- as.data.frame(A$coefficients)
  Ttable$coef.name <- rownames(Ttable)
  Ttable$protein <- transitional_soma$protein[i]
  Ttable$inflection_point <- transitional_soma$Knot[i]
  spline_results_soma2 <- spline_results_soma2 |>
    dplyr::bind_rows(Ttable)

  # prediction and normalization
  # Get estimated marginal means and 95% CI
  # Create a prediction grid
  pred_grid <- expand.grid(
    years_since_pd_onset = seq(-12, 9, by = 0.1),
    group = unique(somalogic_table$group),
    age_at_baseline = mean(somalogic_table$age_at_baseline), # Set to reference or mean value
    gender = "female",
    recent_plasma_donations = "at least one donation",
    cluster = 1
  ) |>
    dplyr::mutate(cluster = as.factor(cluster))

  # Add spline columns
  pred_grid_splines <- lspline::lspline(pred_grid$years_since_pd_onset,
    knots = transitional_soma$Knot[i]
  )
  colnames(pred_grid_splines) <- c("spline1", "spline2")

  pred_grid <- cbind(pred_grid, pred_grid_splines)

  pred_grid$pred <- predict(lmer.linear_spline1, newdata = pred_grid, re.form = NA)

  df_wide <- tidyr::pivot_wider(
    pred_grid |>
      dplyr::select(group, years_since_pd_onset, pred) |>
      dplyr::distinct(),
    names_from = group,
    values_from = pred
  )
  df_wide$normalized <- df_wide$pd - df_wide$control

  df_wide$probe_id <- transitional_soma$protein[i]

  predicted_soma2 <- predicted_soma2 |>
    dplyr::bind_rows(df_wide |>
      dplyr::select(years_since_pd_onset, normalized, probe_id))
}
names(spline_results_soma2)[c(2, 5, 7)] <- c("SE", "p_value", "probe_id")

spline_results_soma2_updated <- spline_results_soma2 |>
  dplyr::filter(coef.name == "lspline::lspline(years_since_pd_onset, transitional_soma$Knot[i])2:grouppd") |>
  dplyr::select(probe_id, Estimate, p_value, inflection_point) |>
  dplyr::inner_join(soma_all |>
    dplyr::select(probe_id, gene_symbol, uniprot_id), by = "probe_id")
predicted_soma2 <- predicted_soma2 |>
  dplyr::inner_join(soma_all |>
    dplyr::select(probe_id, gene_symbol, uniprot_id), by = "probe_id")

### save mixed model and spline mixed model results--results_transitional ####
# mixed model results
all_mixed_model_results2 <- all_combined |>
  dplyr::filter(probe_id %in% transitional_soma$protein)
# head(all_mixed_model_results2)
# linear spline model results
spline_results_combined2 <- spline_results_soma2_updated |>
  dplyr::rename(slope2_group_difference = Estimate, slope2_difference_p_val = p_value)
# combined results
results_transitional <- all_mixed_model_results2 |>
  dplyr::select(platform, probe_id, uniprot_id, gene_symbol, estimate, cohens_d, p_val, q_val) |>
  dplyr::rename(
    log2_FC_all = estimate,
    cohens_d_all = cohens_d,
    p_val_all = p_val,
    q_val_all = q_val
  ) |>
  dplyr::inner_join(spline_results_combined2 |> dplyr::select(-uniprot_id, -gene_symbol), by = "probe_id")


# ==============================================================================
# Group 3a: Progression biomarkers
# ==============================================================================
### Based on linear mixed model results
### Progression biomarkers: identified PD biomarkers with significantly different rate of changes
### from whole sample mixed model results
### Get significant rate of change over all PD phases as in Chronos phase 1

### Somalogic

rate_change_result_v1 <- data.frame()
for (i in 1:length(soma_model)) {
  model_results <- soma_model[[i]]$best_model_summary$tTable |>
    dplyr::filter(coef.names == "years_since_pd_onset:grouppd")
  model_results$protein <- soma_model[[i]]$best_model_summary$protein$ID
  rate_change_result_v1 <- rate_change_result_v1 |>
    dplyr::bind_rows(model_results)
}


rate_change_result_v1 <- rate_change_result_v1 |>
  dplyr::mutate(q_val = p.adjust(p.value, method = "BH")) |>
  dplyr::inner_join(somalogic_tibble$prot_annot |>
    dplyr::select(ID, UniProt, EntrezGeneSymbol), by = c("protein" = "ID"))
rate_change_result_v1 |> dplyr::filter(q_val < 0.05)


progression_soma <- setdiff(
  rate_change_result_v1$protein[rate_change_result_v1$q_val < 0.05],
  c(transitional_soma$protein, transient_preclinical_soma$protein)
)
progression_soma

### Olink
rate_change_result_v2 <- data.frame()
for (i in 1:length(olink_model)) {
  model_results <- olink_model[[i]]$best_model_summary$tTable |>
    dplyr::filter(coef.names == "years_since_pd_onset:grouppd")
  model_results$protein <- olink_model[[i]]$best_model_summary$protein$ID
  rate_change_result_v2 <- rate_change_result_v2 |>
    dplyr::bind_rows(model_results)
}

rate_change_result_v2 <- rate_change_result_v2 |>
  dplyr::mutate(q_val = p.adjust(p.value, method = "BH")) |>
  dplyr::inner_join(olink_tibble$prot_annot |>
    dplyr::select(ID, UNIPROT, ASSAY), by = c("protein" = "ID"))
rate_change_result_v2 |> dplyr::filter(q_val < 0.1) # DDC and PRL
progression_olink <- setdiff(
  rate_change_result_v2$protein[rate_change_result_v2$q_val < 0.05],
  transient_preclinical_olink$protein
)
progression_olink

### Alamar

rate_change_result_v3 <- data.frame()
for (i in 1:length(alamar_model)) {
  model_results <- alamar_model[[i]]$best_model_summary$tTable |>
    dplyr::filter(coef.names == "years_since_pd_onset:grouppd")
  model_results$protein <- alamar_model[[i]]$best_model_summary$protein$ID
  rate_change_result_v3 <- rate_change_result_v3 |>
    dplyr::bind_rows(model_results)
}
# head(alamar_tibble$prot_annot)
rate_change_result_v3 <- rate_change_result_v3 |>
  dplyr::mutate(q_val = p.adjust(p.value, method = "BH")) |>
  dplyr::inner_join(alamar_tibble$prot_annot |>
    dplyr::select(ID, UniProtID, Target), by = c("protein" = "ID"))
rate_change_result_v3 |> dplyr::filter(q_val < 0.05) # DDC
progression_alamar <- setdiff(
  rate_change_result_v3$protein[rate_change_result_v3$q_val < 0.05],
  transient_preclinical_alamar$protein
)

# ------------------------------------------------------------------------------
# Save model results for progression biomarkers
# ------------------------------------------------------------------------------
progression_probe_id <- c(progression_soma, progression_olink)
progression_biomarker_table <- rate_change_result_v1 |>
  dplyr::filter(protein %in% progression_probe_id) |>
  dplyr::rename(
    probe_id = protein,
    uniprot_id = UniProt,
    gene_symbol = EntrezGeneSymbol
  ) |>
  dplyr::bind_rows(rate_change_result_v2 |>
    dplyr::filter(protein %in% progression_probe_id) |>
    dplyr::rename(
      probe_id = protein,
      uniprot_id = UNIPROT,
      gene_symbol = ASSAY
    )) |>
  dplyr::select(probe_id, uniprot_id, gene_symbol, Estimate, cohensD, p.value, q_val) |>
  dplyr::rename(
    slope_group_difference = Estimate,
    cohens_D = cohensD,
    p_val = p.value
  )

progression_biomarker_table


# ==============================================================================
# Group 3b: Prodromal biomarkers
# ==============================================================================
### Summary statistics from linear mixed model and joint model
### Prodromal biomarkers: significantly different between PD and control in whole sample linear mixed model
### from at least one platform and were early PD biomarkers based on preDx samples linear mixed models
### or candidate risk PD biomarkers based on joint modeling, and they do not belong to above 3 groups

prodromal_biomarker_soma <- setdiff(
  intersect(
    soma_all$probe_id[soma_all$q_val < 0.05],
    union(
      soma_pre$probe_id[soma_pre$q_val < 0.05],
      soma_jm$probe_id[soma_jm$p_val < 0.01 & soma_jm$rhat < 1.05]
    )
  ),
  c(transient_preclinical_soma$protein, transitional_soma$protein, progression_soma)
)

prodromal_biomarker_olink <- setdiff(
  unique(intersect(
    olink_all$probe_id[olink_all$q_val < 0.05],
    union(
      olink_pre$probe_id[olink_pre$q_val < 0.05],
      olink_jm$probe_id[olink_jm$p_val < 0.01 & olink_jm$rhat < 1.05]
    )
  )),
  c(transient_preclinical_olink$protein, progression_olink)
)


prodromal_biomarker_alamar <- setdiff(
  unique(intersect(
    alamar_all$probe_id[alamar_all$q_val < 0.05],
    union(
      alamar_pre$probe_id[alamar_pre$q_val < 0.05],
      alamar_jm$probe_id[alamar_jm$p_val < 0.01 & alamar_jm$rhat < 1.05]
    )
  )),
  c(transient_preclinical_alamar$protein)
)
prodromal_biomarker_alamar


prodromal_probe_id <- setdiff(
  c(prodromal_biomarker_soma, prodromal_biomarker_olink, prodromal_biomarker_alamar),
  c(transient_preclinical_probe_id, transitional_probe_id, progression_probe_id)
)
prodromal_probe_id

# ------------------------------------------------------------------------------
# Save mixed model and joint model results for prodromal biomarkers
# ------------------------------------------------------------------------------

results_prodromal <- all_combined |>
  dplyr::filter(probe_id %in% prodromal_probe_id) |>
  dplyr::rename(
    log2_FC_all = estimate,
    cohens_d_all = cohens_d,
    p_val_all = p_val,
    q_val_all = q_val
  ) |>
  dplyr::inner_join(
    preDx_combined |>
      dplyr::select(-uniprot_id, -gene_symbol, -platform) |>
      dplyr::rename(
        log2_FC_pre = estimate,
        cohens_d_pre = cohens_d,
        p_val_pre = p_val,
        q_val_pre = q_val
      ),
    by = "probe_id"
  ) |>
  dplyr::left_join(
    risk_combined |>
      dplyr::select(probe_id, risk, risk_lower, risk_upper, p_val) |>
      dplyr::rename(p_val_risk = p_val),
    by = "probe_id"
  )

# remove probe ids where their uniprotID was in transient preclinical, transitional or progression biomarkers
remove_protein_list <- c(
  results_transient_preclinical$uniprot_id,
  results_transitional$uniprot_id,
  progression_biomarker_table$uniprot_id
)
length(unique(remove_protein_list))
results_prodromal |>
  dplyr::filter(uniprot_id %in% remove_protein_list) # only DDC
results_prodromal <- results_prodromal |> dplyr::filter(uniprot_id %nin% remove_protein_list)


# ==============================================================================
# Group 3c: Post-diagnosis biomarkers
# ==============================================================================
### The PD signals that were not in either of the above four groups
### but significantly different between PD and control in postDx samples or significantly different
### in the second slope with inflection point after PD onset from linear spline mixed model results

#### any significant different rate of change from the second slope (non-linear models)
postDx_soma_ls <- somalogic_all_ls_selected_BIC |>
  dplyr::filter(Knot > 0 & group_diff_sig %in% c("second", "both"))
postDx_soma_ls # 2
postDx_olink_ls <- olink_all_ls_selected_BIC |>
  dplyr::filter(Knot > 0 & group_diff_sig %in% c("second", "both"))
postDx_olink_ls # none
postDx_alamar_ls <- alamar_all_ls_selected_BIC |>
  dplyr::filter(Knot > 0 & group_diff_sig %in% c("second", "both"))
postDx_alamar_ls # none

# somalogic results from mixed models, joint models
somalogic_results <- soma_post |>
  dplyr::inner_join(soma_pre |>
    dplyr::select(p_val, probe_id) |>
    dplyr::rename(p_val_pre = p_val), by = "probe_id") |>
  dplyr::inner_join(soma_jm |>
    dplyr::select(p_val, probe_id) |>
    dplyr::rename(p_val_jm = p_val), by = "probe_id") |>
  dplyr::left_join(postDx_soma_ls |>
    dplyr::select(protein, Knot) |>
    dplyr::rename(probe_id = protein), by = "probe_id")
somalogic_results |>
  dplyr::filter(!is.na(Knot))

postDx_biomarkers_soma <- somalogic_results |>
  dplyr::filter((abs(cohens_d) > 0.2 & p_val < 0.005 | !is.na(Knot)) & p_val_pre > 0.05 & p_val_jm > 0.05 &
    probe_id %nin% c(transitional_probe_id, transient_preclinical_probe_id, progression_probe_id, prodromal_probe_id)) |>
  dplyr::select(probe_id, uniprot_id, gene_symbol, estimate, cohens_d, p_val, q_val, p_val_pre, p_val_jm, Knot) |>
  dplyr::arrange(p_val)
postDx_biomarkers_soma |> dplyr::filter(q_val >= 0.05)
# olink results from mixed models, joint models
olink_results <- olink_post |>
  dplyr::inner_join(olink_pre |>
    dplyr::select(p_val, probe_id) |>
    dplyr::rename(p_val_pre = p_val), by = "probe_id") |>
  dplyr::inner_join(olink_jm |>
    dplyr::select(p_val, probe_id) |>
    dplyr::rename(p_val_jm = p_val), by = "probe_id")

postDx_biomarkers_olink <- olink_results |>
  dplyr::filter(abs(cohens_d) > 0.2 & p_val < 0.005 & p_val_pre > 0.05 & p_val_jm > 0.05 &
    probe_id %nin% c(transitional_probe_id, transient_preclinical_probe_id, progression_probe_id, prodromal_probe_id)) |>
  dplyr::select(probe_id, uniprot_id, gene_symbol, estimate, cohens_d, p_val, q_val, p_val_pre, p_val_jm) |>
  dplyr::arrange(p_val)


# olink results from mixed models, joint models
alamar_results <- alamar_post |>
  dplyr::inner_join(alamar_pre |>
    dplyr::select(p_val, probe_id) |>
    dplyr::rename(p_val_pre = p_val), by = "probe_id") |>
  dplyr::inner_join(alamar_jm |>
    dplyr::select(p_val, probe_id) |>
    dplyr::rename(p_val_jm = p_val), by = "probe_id") |>
  dplyr::left_join(postDx_alamar_ls |>
    dplyr::select(protein, Knot) |>
    dplyr::rename(probe_id = protein), by = "probe_id")
# head(alamar_results)
postDx_biomarkers_alamar <- alamar_results |>
  dplyr::filter((abs(cohens_d) > 0.2 & p_val < 0.005 | !is.na(Knot)) & p_val_pre > 0.05 & p_val_jm > 0.05 &
    probe_id %nin% c(transitional_probe_id, transient_preclinical_probe_id, progression_probe_id, prodromal_probe_id)) |>
  dplyr::select(probe_id, uniprot_id, gene_symbol, estimate, cohens_d, p_val, q_val, p_val_pre, p_val_jm) |>
  dplyr::arrange(p_val)

# ------------------------------------------------------------------------------
# Save overall model results for post-diagnosis biomarkers
# ------------------------------------------------------------------------------

postDx_PD_probe_id <- c(postDx_biomarkers_soma$probe_id, postDx_biomarkers_olink$probe_id)
postDx_PD_probe_id
postDx_biomarkers_table <- postDx_biomarkers_soma |>
  dplyr::mutate(platform = "Somalogic11k") |>
  dplyr::bind_rows(postDx_biomarkers_olink |>
    dplyr::mutate(platform = "OlinkHT"))


# ==============================================================================
# Save all predicted values for visualization
# ==============================================================================
# Save all predictions based on linear spline mixed models (transient preclinical and transitional biomarkers)
predicted_results_1 <- predicted_soma |>
  dplyr::mutate(
    platform = "SomaScanV5",
    biomarker_type = "transient_preclinical"
  ) |>
  dplyr::bind_rows(predicted_soma2 |>
    dplyr::mutate(
      platform = "SomaScanV5",
      biomarker_type = "transitional"
    )) |>
  dplyr::bind_rows(predicted_olink |>
    dplyr::mutate(
      platform = "Olink_HT",
      biomarker_type = "transient_preclinical"
    )) |>
  dplyr::bind_rows(predicted_alamar |>
    dplyr::mutate(
      platform = "Alamarbio_CNS",
      biomarker_type = "transient_preclinical"
    ))

# ------------------------------------------------------------------------------
# Create and save predictions based on linear mixed models
# (progression, prodromal and postDx biomarkers)
# ------------------------------------------------------------------------------
#### Prediction and normalization
ctrl <- nlme::lmeControl(opt = "optim")

predicted_results <- data.frame()
# somalogic
for (i in c(prodromal_biomarker_soma, progression_soma, postDx_biomarkers_soma$probe_id)) {
  # print(i)
  mydataset <- somalogic_table |>
    dplyr::select(all_of(i), years_since_pd_onset, group, age_at_baseline, gender, recent_plasma_donations, cluster, patient_id)
  names(mydataset)[1] <- "protein"
  model <- nlme::lme(
    fixed = protein ~ years_since_pd_onset * group + age_at_baseline + gender +
      recent_plasma_donations + cluster,
    random = ~ 1 + years_since_pd_onset | patient_id,
    data = mydataset,
    control = ctrl
  )

  # Get estimated marginal means and 95% CI
  em <- emmeans::emmeans(model, ~ years_since_pd_onset * group, at = list(years_since_pd_onset = seq(-12, 9, by = 0.1)))

  # Convert emmeans object to a data frame
  pred_df <- as.data.frame(em)
  # head(pred_df)

  df_wide <- tidyr::pivot_wider(
    pred_df |>
      dplyr::select(group, years_since_pd_onset, emmean, SE),
    names_from = group,
    values_from = c(emmean, SE)
  )

  # using subtraction
  df_wide$normalized <- df_wide$emmean_pd - df_wide$emmean_control
  df_wide$SE_normalized <- sqrt(
    (df_wide$SE_pd)^2 +
      (df_wide$SE_control)^2
  )

  df_wide$lower_CI <- df_wide$normalized - 1.96 * df_wide$SE_normalized
  df_wide$upper_CI <- df_wide$normalized + 1.96 * df_wide$SE_normalized

  df_wide$probe_id <- i
  df_wide$platform <- "SomaScanV5"

  predicted_results <- predicted_results |>
    dplyr::bind_rows(df_wide)
}

# olink
for (i in c(prodromal_biomarker_olink, progression_olink, postDx_biomarkers_olink$probe_id)) {
  # print(i)
  mydataset <- olink_table |>
    dplyr::select(all_of(i), years_since_pd_onset, group, age_at_baseline, gender, recent_plasma_donations, cluster, patient_id)
  names(mydataset)[1] <- "protein"
  model <- nlme::lme(
    fixed = protein ~ years_since_pd_onset * group + age_at_baseline + gender +
      recent_plasma_donations + cluster,
    random = ~ 1 + years_since_pd_onset | patient_id,
    data = mydataset,
    control = ctrl
  )

  # Get estimated marginal means and 95% CI
  em <- emmeans::emmeans(model, ~ years_since_pd_onset * group, at = list(years_since_pd_onset = seq(-12, 9, by = 0.1)))

  # Convert emmeans object to a data frame
  pred_df <- as.data.frame(em)
  # head(pred_df)

  df_wide <- tidyr::pivot_wider(
    pred_df |>
      dplyr::select(group, years_since_pd_onset, emmean, SE),
    names_from = group,
    values_from = c(emmean, SE)
  )

  # using subtraction--this one
  df_wide$normalized <- df_wide$emmean_pd - df_wide$emmean_control
  df_wide$SE_normalized <- sqrt(
    (df_wide$SE_pd)^2 +
      (df_wide$SE_control)^2
  )

  df_wide$lower_CI <- df_wide$normalized - 1.96 * df_wide$SE_normalized
  df_wide$upper_CI <- df_wide$normalized + 1.96 * df_wide$SE_normalized

  df_wide$probe_id <- i
  df_wide$platform <- "Olink_HT"

  predicted_results <- predicted_results |>
    dplyr::bind_rows(df_wide)
}

# alamarbio
for (i in c(setdiff(prodromal_biomarker_alamar, "DDC"), postDx_biomarkers_alamar$probe_id)) {
  # print(i)
  mydataset <- alamar_table |>
    dplyr::select(all_of(i), years_since_pd_onset, group, age_at_baseline, gender, recent_plasma_donations, cluster, patient_id)
  names(mydataset)[1] <- "protein"
  model <- nlme::lme(
    fixed = protein ~ years_since_pd_onset * group + age_at_baseline + gender +
      recent_plasma_donations + cluster,
    random = ~ 1 + years_since_pd_onset | patient_id,
    data = mydataset,
    control = ctrl
  )

  # Get estimated marginal means and 95% CI
  em <- emmeans::emmeans(model, ~ years_since_pd_onset * group, at = list(years_since_pd_onset = seq(-12, 9, by = 0.1)))

  # Convert emmeans object to a data frame
  pred_df <- as.data.frame(em)
  # head(pred_df)

  df_wide <- tidyr::pivot_wider(
    pred_df |>
      dplyr::select(group, years_since_pd_onset, emmean, SE),
    names_from = group,
    values_from = c(emmean, SE)
  )

  # using subtraction--this one
  df_wide$normalized <- df_wide$emmean_pd - df_wide$emmean_control
  df_wide$SE_normalized <- sqrt(
    (df_wide$SE_pd)^2 +
      (df_wide$SE_control)^2
  )

  df_wide$lower_CI <- df_wide$normalized - 1.96 * df_wide$SE_normalized
  df_wide$upper_CI <- df_wide$normalized + 1.96 * df_wide$SE_normalized

  df_wide$probe_id <- i
  df_wide$platform <- "Alamarbio_CNS"

  predicted_results <- predicted_results |>
    dplyr::bind_rows(df_wide)
}

predicted_results_2 <- predicted_results |>
  dplyr::mutate(biomarker_type = ifelse(
    probe_id %in% progression_probe_id, "progression",
    ifelse(probe_id %in% prodromal_probe_id, "prodromal", "postDx")
  )) |>
  dplyr::inner_join(all_combined |>
    dplyr::select(probe_id, gene_symbol, uniprot_id), by = "probe_id")

# ------------------------------------------------------------------------------
# Combine all normalized mean predictions
# ------------------------------------------------------------------------------
combined_predicted_results <- predicted_results_1 |>
  dplyr::select(-platform) |>
  dplyr::bind_rows(predicted_results_2 |>
    dplyr::select(years_since_pd_onset, normalized, probe_id, gene_symbol, uniprot_id, biomarker_type))
## save biomarker groups ####
biomarker_groups <- combined_predicted_results |>
  dplyr::select(probe_id, gene_symbol, uniprot_id, biomarker_type) |>
  dplyr::distinct()


# ==============================================================================
# Examine medication effects on PD biomarker groups
# ==============================================================================

Olink_med_on_biomarkers <- Olink_med_table_path |>
  arrow::open_dataset() |>
  dplyr::collect() |>
  dplyr::left_join(biomarker_groups, by = c("protein_id" = "probe_id")) |>
  dplyr::rename(gene_symbol = gene_symbol.x) |>
  dplyr::filter(!is.na(biomarker_type))
olink_med_list <- Olink_med_on_biomarkers |> dplyr::filter(abs(cohensD) >= 0.8 & p.value < 0.05)
olink_med_list |> dplyr::select(protein_id, gene_symbol, cohensD, estimate, p.value) # DDC


soma_med_on_biomarkers <- SomaLogic_med_table_path |>
  arrow::open_dataset() |>
  dplyr::collect() |>
  dplyr::left_join(biomarker_groups, by = c("protein_id" = "probe_id")) |>
  dplyr::rename(gene_symbol = gene_symbol.x) |>
  dplyr::filter(!is.na(biomarker_type))
# head(soma_med_on_biomarkers)
soma_med_list <- soma_med_on_biomarkers |> dplyr::filter(abs(cohensD) >= 0.8 & p.value < 0.05)
soma_med_list |> dplyr::select(protein_id, gene_symbol, cohensD, estimate, p.value)


Alamar_med_on_biomarkers <- Alamar_med_table_path |>
  arrow::open_dataset() |>
  dplyr::collect() |>
  dplyr::left_join(biomarker_groups, by = c("protein_id" = "probe_id")) |>
  dplyr::rename(gene_symbol = gene_symbol.x) |>
  dplyr::filter(!is.na(biomarker_type))
# head(soma_med_on_biomarkers)
Alamar_med_on_biomarkers |> dplyr::filter(abs(cohensD) >= 0.8 & p.value < 0.05) #

medication_list_gene_symbol <- c(olink_med_list$gene_symbol, soma_med_list$gene_symbol)


# ==============================================================================
# Remove medication effects from PD biomarker groups
# ==============================================================================
progression_list <- progression_biomarker_table |> dplyr::filter(gene_symbol %nin% medication_list_gene_symbol)
transitional_list <- results_transitional |> dplyr::filter(gene_symbol %nin% medication_list_gene_symbol)
postDx_list <- postDx_biomarkers_table |> dplyr::filter(gene_symbol %nin% medication_list_gene_symbol)

transient_preclinical_list <- results_transient_preclinical |> dplyr::filter(gene_symbol %nin% medication_list_gene_symbol)
prodromal_list <- results_prodromal |> dplyr::filter(gene_symbol %nin% medication_list_gene_symbol)


# ==============================================================================
# Create final 3-group PD trajectory tables for manuscript
# ==============================================================================
# Group 1: Transient preclinical
group1_early_nonlinear <- transient_preclinical_list

# Group 2: Prodromal
group2_prodromal_early_linear <- prodromal_list |>
  dplyr::select(
    platform, probe_id, uniprot_id, gene_symbol,
    log2_FC_all, cohens_d_all, p_val_all, q_val_all,
    log2_FC_pre, cohens_d_pre, p_val_pre, q_val_pre,
    risk, risk_lower, risk_upper, p_val_risk
  )


# Group 3: Late and progression
# Save postDx results; for progression and transitional, add rate of change results
group3_transitional <- soma_post |>
  dplyr::filter(probe_id %in% transitional_list$probe_id) |>
  dplyr::mutate(platform = "SomaScanV5") |>
  dplyr::rename(
    log2_FC_post = estimate,
    cohens_d_post = cohens_d,
    p_val_post = p_val,
    q_val_post = q_val
  ) |>
  dplyr::select(platform, probe_id, uniprot_id, gene_symbol, log2_FC_post, cohens_d_post, p_val_post, q_val_post) |>
  dplyr::inner_join(transitional_list |>
    dplyr::select(probe_id, slope2_group_difference, slope2_difference_p_val, inflection_point) |>
    dplyr::rename(
      slope_group_difference = slope2_group_difference,
      slope_difference_p_val = slope2_difference_p_val
    ), by = "probe_id")


group3_progression <- soma_post |>
  dplyr::filter(probe_id %in% progression_list$probe_id) |>
  dplyr::mutate(platform = "SomaScanV5") |>
  dplyr::rename(
    log2_FC_post = estimate,
    cohens_d_post = cohens_d,
    p_val_post = p_val,
    q_val_post = q_val
  ) |>
  dplyr::select(platform, probe_id, uniprot_id, gene_symbol, log2_FC_post, cohens_d_post, p_val_post, q_val_post) |>
  dplyr::bind_rows(olink_post |>
    dplyr::filter(probe_id %in% progression_list$probe_id) |>
    dplyr::mutate(platform = "Olink_HT") |>
    dplyr::rename(
      log2_FC_post = estimate,
      cohens_d_post = cohens_d,
      p_val_post = p_val,
      q_val_post = q_val
    ) |>
    dplyr::select(platform, probe_id, uniprot_id, gene_symbol, log2_FC_post, cohens_d_post, p_val_post, q_val_post)) |>
  dplyr::inner_join(progression_list |>
    dplyr::select(probe_id, slope_group_difference, p_val) |>
    dplyr::rename(slope_difference_p_val = p_val), by = "probe_id")

group3_postDx <- soma_post |>
  dplyr::filter(probe_id %in% postDx_list$probe_id) |>
  dplyr::mutate(platform = "SomaScanV5") |>
  dplyr::rename(
    log2_FC_post = estimate,
    cohens_d_post = cohens_d,
    p_val_post = p_val,
    q_val_post = q_val
  ) |>
  dplyr::select(platform, probe_id, uniprot_id, gene_symbol, log2_FC_post, cohens_d_post, p_val_post, q_val_post) |>
  dplyr::bind_rows(olink_post |>
    dplyr::filter(probe_id %in% postDx_list$probe_id) |>
    dplyr::mutate(platform = "Olink_HT") |>
    dplyr::rename(
      log2_FC_post = estimate,
      cohens_d_post = cohens_d,
      p_val_post = p_val,
      q_val_post = q_val
    ) |>
    dplyr::select(
      platform, probe_id, uniprot_id, gene_symbol,
      log2_FC_post, cohens_d_post, p_val_post, q_val_post
    ))


group3_late_progression <- group3_transitional |>
  dplyr::bind_rows(group3_progression) |>
  dplyr::bind_rows(group3_postDx)


# ==============================================================================
# Generate Supplementary Table 7
# ==============================================================================

ST7 <- group1_early_nonlinear |>
  dplyr::mutate(
    "Trajectory group" = "inflection years before onset"
  ) |>
  dplyr::bind_rows(
    group2_prodromal_early_linear |>
      dplyr::mutate(
        "Trajectory group" = "early linear change"
      )
  ) |>
  dplyr::bind_rows(
    group3_late_progression |>
      dplyr::mutate(
        "Trajectory group" = "late linear and nonlinear change"
      )
  )

ST7 <- dplyr::bind_rows(
  group1_early_nonlinear |>
    dplyr::mutate(
      biomarker_category = "inflection years before onset"
    ),
  group2_prodromal_early_linear |>
    dplyr::mutate(
      biomarker_category = "early linear change"
    ),
  group3_late_progression |>
    dplyr::mutate(
      biomarker_category = "late linear and nonlinear change"
    )
) |>
  dplyr::select(
    "Trajectory group" = "biomarker_category",
    "Proteomics Platform" = "platform",
    "Probe ID" = "probe_id",
    "UniProt ID" = "uniprot_id",
    "Gene Symbol" = "gene_symbol",
    "Log2 FC linear mixed model all samples" = "log2_FC_all",
    "Cohen's d all samples" = "cohens_d_all",
    "Pvalue linear mixed model all samples" = "p_val_all",
    "Qvalue linear mixed model all samples" = "q_val_all",
    "Difference in the rate of change between PD and control at the first slope linear spline mixed model" = "slope1_group_difference",
    "Pvalue at the first slope linear spline mixed model" = "slope1_difference_p_val",
    "Inflection point linear spline mixed model" = "inflection_point",
    "Log2 FC linear mixed model pre-onset samples" = "log2_FC_pre",
    "Cohen's d pre-onset samples" = "cohens_d_pre",
    "Pvalue linear mixed model pre-onset samples" = "p_val_pre",
    "Qvalue linear mixed model pre-onset samples" = "q_val_pre",
    "Mean risk per 0.1 unit change joint modeling" = "risk",
    "Lower limit of risk joint modeling" = "risk_lower",
    "Upper limit of risk joint modeling" = "risk_upper",
    "Pvalue joint modeling" = "p_val_risk",
    "Log2 FC linear mixed model post-onset samples" = "log2_FC_post",
    "Cohen's d post-onset samples" = "cohens_d_post",
    "Pvalue linear mixed model post-onset samples" = "p_val_post",
    "Qvalue linear mixed model post-onset samples" = "q_val_post",
    "Slope difference linear mixed model all samples" = "slope_group_difference",
    "Pvalue slope difference linear mixed model all samples" = "slope_difference_p_val"
  )

ST7 |>
  readr::write_tsv(output_file)
