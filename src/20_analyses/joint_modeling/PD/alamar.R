## SPDX-FileCopyrightText: 2026 Alkahest, Inc. <chronospd_data_access@alkahest.com>
##
## SPDX-License-Identifier: MIT

here::i_am("src/20_analyses/joint_modeling/PD/alamar.R")

##### !!!!! #####
# Run 1
# To use slurm: convert the integer to take arguments
protein_index <- 1

future::plan(future::multisession, workers = parallelly::availableCores())

# screen model ####
input_phase1_critical_variables_file <- glue::glue("{get_config()$input_dir}/name=Chronos_phase1_critical_variables_table/version=29432/dataset.parquet")
input_phase1_RWD_summary_file <- glue::glue("{get_config()$input_dir}/name=phase1_claim_summary/version=none/dataset.parquet")
input_alamar_results_file <- glue::glue("{get_config()$input_dir}/name=chronos_alamar_merged-pd_master_summary_2/version=29830/dataset.parquet")
input_alamar_tibble_file <- glue::glue("{get_config()$input_dir}/name=Chronos_phase1_proteomics_alamar_analysis_ready_three_tibble/version=29464/dataset.rds")

output_result_file <- glue::glue("{get_config()$checkpoint_dir}/joint_modeling/PD/alamar.screening.protein.{protein_index}.rds")
dir.create(dirname(output_result_file), showWarnings = FALSE, recursive = TRUE)

"%nin%" <- Negate("%in%")
### general data loading ####
# critical variable table ####
phase1_critical_variables <- input_phase1_critical_variables_file |>
  arrow::open_dataset() |>
  dplyr::collect()

# RWD summary data from claim ####
results_phase1 <- input_phase1_RWD_summary_file |>
  arrow::open_dataset() |>
  dplyr::select(patient_id, first_date_id, last_date_id) |>
  dplyr::distinct() |>
  dplyr::collect()

# sample-level data with PD onset time ####
sample_PD_onset <- results_phase1 |>
  dplyr::left_join(phase1_critical_variables, by = "patient_id") |>
  dplyr::mutate(
    Event_date = ifelse(group == "pd", estimated_pd_onset_date, last_date_id),
    min_sample_date = min(sample_collection_date),
    birth_year = as.Date(paste0(birth_year, "-01-01"))
  ) |>
  dplyr::filter(sample_collection_date < Event_date) |>
  dplyr::group_by(patient_id) |>
  dplyr::mutate(
    n_sample = dplyr::n()
  ) |>
  dplyr::ungroup() |>
  dplyr::mutate(
    time_protein = (as.numeric(difftime(as.Date(sample_collection_date), # years since "baseline" of the samples
                                        as.Date(min_sample_date),
                                        units = "days"
    )) + 1) / 365.25,
    time_event = as.numeric(difftime(as.Date(Event_date, origin = "1970-01-01"), # years since "baseline of the event
                                     as.Date(min_sample_date),
                                     units = "days"
    )) / 365.25,
    age_event = as.numeric(difftime(as.Date(Event_date, origin = "1970-01-01"), # age at event
                                    as.Date(birth_year),
                                    units = "days"
    )) / 365.25,
    age_entry = as.numeric(difftime(as.Date(first_date_id, origin = "1970-01-01"), # age with the first medical record (claim and medication)
                                    as.Date(birth_year),
                                    units = "days"
    )) / 365.25,
    age_protein = (as.numeric(difftime(as.Date(sample_collection_date), # age at donation, so age for the proteins measured
                                       as.Date(birth_year),
                                       units = "days"
    )) + 1) / 365.25,
  )

### alamar JM ####
# alamar data ####
# model results:
alamar_results <- input_alamar_results_file |>
  arrow::open_dataset() |>
  dplyr::collect()

alamar_detection <- alamar_results |>
  dplyr::filter(q_val_all < 0.05) # 14

# proteome:
alamar_tibble <- readRDS(input_alamar_tibble_file)
alamar_table <- alamar_tibble$meta_data |>
  dplyr::inner_join(
    alamar_tibble$prot_data,
    by = "ID"
  )

# data preparing ####
# merge sample-level data with protein read:
sample_protein_PD_onset_alamar <- sample_PD_onset |>
  dplyr::inner_join(
    alamar_table |>
      dplyr::select(PARENT_BARCODE, tidyselect::all_of(alamar_results$probe_id)),
    by = c("sample_id" = "PARENT_BARCODE")
  )

# survival data:
survival_data_alamar <- sample_protein_PD_onset_alamar |>
  dplyr::select(
    patient_id, time_event, age_event, age_entry, Event_date, group, gender,
    race, birth_year
  ) |>
  dplyr::distinct() |>
  dplyr::mutate(
    Time = time_event,
    Time_age = age_event,
    Status = ifelse(group == "pd", 1, 0)
  )

# JM ####
# data list:
data_list_alamar <- list(
  longitudinal = sample_protein_PD_onset_alamar,
  survival = survival_data_alamar
)

# JM specific, using age as time and add entry time, not paired:
jmdata_alamar <- joint_data(
  longitudinal = data_list_alamar$longitudinal |>
    dplyr::filter(age_entry < age_event),
  survival = data_list_alamar$survival |>
    dplyr::filter(age_entry < age_event),
  id_var = patient_id,
  time_var = age_protein
)

spec <- JointModelSpec$new(
  longitudinal_formula = ~ gender + cluster + recent_plasma_donations + age_protein,
  survival_formula = survival::Surv(age_entry, Time_age, Status) ~ gender,
  random_formula = ~ 1 + age_protein | patient_id,
  association = "value"
)

spec$set_mcmc(
  chains = 4,
  iter = 30000,
  warmup = 3000,
  thin = 2
)

protein_subset <- unique(alamar_detection$probe_id)[protein_index]

results <- fit_joint_models(
  data = jmdata_alamar,
  spec = spec,
  proteins = protein_subset,
  progress = TRUE
)

results |>
  readr::write_rds(output_result_file)

future::plan(future::sequential)
