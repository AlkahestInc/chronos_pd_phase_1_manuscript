## SPDX-FileCopyrightText: 2026 Alkahest, Inc. <chronospd_data_access@alkahest.com>
##
## SPDX-License-Identifier: MIT

here::i_am("src/20_analyses/joint_modeling/comorbidity/model_processing/Somalogic_PD_comorbidity.R")

# JM functions are loaded automatically from packages/manuscript via .Rprofile
# data loading ####
cli::cli_h1("Somalogic PD Comorbidity Joint Modeling")
cli::cli_alert_info("Loading input data...")

protein_list <- arrow::read_parquet(glue::glue("{get_config()$checkpoint_dir}/name=protein_list/version=none/dataset.parquet"))
JM_lme_somalogic <- arrow::read_parquet(glue::glue("{get_config()$checkpoint_dir}/name=JM_lme_somalogic/version=none/dataset.parquet"))
JM_surv_somalogic <- arrow::read_parquet(glue::glue("{get_config()$checkpoint_dir}/name=JM_surv_somalogic/version=none/dataset.parquet"))

# Get all comorbidity codes (3-digit ICD10 codes) in this dataset
diagnosis_code_list <- sort(unique(JM_surv_somalogic$diagnosis_code_3))
n_diag <- length(diagnosis_code_list)

# SomaScan proteins only
soma_idx      <- which(protein_list$platform == "SomaScanV5")
soma_proteins <- protein_list$probe_id[soma_idx]
n_prot        <- length(soma_proteins)

cli::cli_alert_info(glue::glue(
  "Found {n_diag} diagnosis codes and {n_prot} SomaScan proteins"
))
cli::cli_alert_info(glue::glue(
  "Total combinations to process: {n_diag * n_prot}"
))

# Setup output directory
output_base_full <- glue::glue("{get_config()$checkpoint_dir}/joint_modeling/comorbidity/pd_comorbidity/full_results")
output_base_summary <- glue::glue("{get_config()$checkpoint_dir}/joint_modeling/comorbidity/pd_comorbidity/summary_only")
dir.create(output_base_full, showWarnings = FALSE, recursive = TRUE)
dir.create(output_base_summary, showWarnings = FALSE, recursive = TRUE)

# JM model specification (same for all combinations)
spec <- JointModelSpec$new(
  longitudinal_formula = ~ age_protein + gender + cluster + recent_plasma_donations,
  survival_formula = survival::Surv(age_entry, Time_age, Status) ~ gender,
  random_formula = ~ 1 + age_protein | PATIENT_ID,
  association = "value"
)
spec$set_mcmc(
  chains = 4,
  iter = 30000,
  warmup = 3000,
  thin = 2
)

options(future.rng.onMisuse = "ignore")
future::plan(
  future::multisession,
  workers = future::availableCores() - 1
)

# Loop through all diagnosis × protein combinations ####
cli::cli_h2("Starting joint modeling analysis")

completed <- 0
skipped <- 0
failed <- 0
total <- n_diag * n_prot

for (diag_index in seq_along(diagnosis_code_list)) {
  diag_code <- diagnosis_code_list[diag_index]

  for (protein_index in seq_along(soma_proteins)) {
    protein_subset <- soma_proteins[protein_index]

    # Define output filepaths
    output_file_full <- glue::glue("{output_base_full}/somalogic_JM.{diag_code}.{protein_index}.{protein_subset}.RDS")
    output_file_summary <- glue::glue("{output_base_summary}/somalogic_JM.{diag_code}.{protein_index}.{protein_subset}.RDS")

    # Skip if already completed
    if (file.exists(output_file_full) && file.exists(output_file_summary)) {
      skipped <- skipped + 1
      cli::cli_alert_info(glue::glue(
        "[{completed + skipped + failed}/{total}] Skipping {diag_code} × {protein_subset} (already exists)"
      ))
      next
    }

    cli::cli_alert_info(glue::glue(
      "[{completed + skipped + failed + 1}/{total}] Running {diag_code} × {protein_subset}"
    ))

    # Try to run the model, catch errors to continue if one fails
    tryCatch({
      # Filter to this diagnosis and report sizes
      long_sub <- JM_lme_somalogic  |> dplyr::filter(diagnosis_code_3 == diag_code, age_entry < age_event)
      surv_sub <- JM_surv_somalogic |> dplyr::filter(diagnosis_code_3 == diag_code, age_entry < age_event)

      if (nrow(long_sub) == 0L || nrow(surv_sub) == 0L) {
        cli::cli_alert_warning(glue::glue(
          "No data after filtering for {diag_code}; skipping"
        ))
        skipped <- skipped + 1
        next
      }

      # Create joint_data
      data_list_somalogic <- list(
        longitudinal = long_sub,
        survival     = surv_sub
      )
      jmdata_somalogic <- joint_data(
        longitudinal = data_list_somalogic$longitudinal,
        survival     = data_list_somalogic$survival,
        id_var   = PATIENT_ID,
        time_var = age_protein
      )

      # Fit joint model
      results <- fit_joint_models(
        jmdata_somalogic,
        spec,
        proteins = protein_subset,
        parallel = FALSE,
        progress = FALSE
      )

      # Save results
      saveRDS(results, output_file_full)
      saveRDS(results$summaries, output_file_summary)

      completed <- completed + 1
      cli::cli_alert_success(glue::glue(
        "Completed {diag_code} × {protein_subset}"
      ))

    }, error = function(e) {
      failed <- failed + 1
      cli::cli_alert_danger(glue::glue(
        "Failed {diag_code} × {protein_subset}: {e$message}"
      ))
    })
  }
}

# Final summary
cli::cli_h2("Analysis Complete")
cli::cli_alert_success(glue::glue("Completed: {completed}"))
cli::cli_alert_info(glue::glue("Skipped: {skipped}"))
if (failed > 0) {
  cli::cli_alert_warning(glue::glue("Failed: {failed}"))
}
