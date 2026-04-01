## SPDX-FileCopyrightText: 2026 Alkahest, Inc. <chronospd_data_access@alkahest.com>
##
## SPDX-License-Identifier: MIT

here::i_am("src/20_analyses/joint_modeling/comorbidity/model_processing/joint_modeling_preparation.R")
"%nin%" <- Negate("%in%")

# Setup output directory ####
# Note: These files should eventually be published to the data publishing repository
# and made available via get_config()$input_dir for the modeling scripts
data_prep_base <- glue::glue("{get_config()$checkpoint_dir}/joint_modeling/comorbidity/model_processing")
dir.create(data_prep_base, showWarnings = FALSE, recursive = TRUE)

# Helper function to save datasets in versioned parquet format
save_versioned_dataset <- function(data, dataset_name, base_dir = data_prep_base, version = "none") {
  output_path <- glue::glue("{base_dir}/name={dataset_name}/version={version}")
  dir.create(output_path, showWarnings = FALSE, recursive = TRUE)
  arrow::write_parquet(data, glue::glue("{output_path}/dataset.parquet"))
  cli::cli_alert_success(glue::glue("Saved: {dataset_name}"))
}

cli::cli_alert_info(glue::glue("Output directory: {data_prep_base}"))

#####
# critical variable ####
#####
phase1_critical_variables <- glue::glue(
  "{get_config()$input_dir}/name={name}/version={version}/dataset.parquet",
  name = "Chronos_phase1_critical_variables_table",
  version = 29432
) |>
  arrow::open_dataset() |>
  dplyr::collect()

#####
# somalogic data ####
#####
somalogic_tibble <- glue::glue(
  "{get_config()$input_dir}/name={name}/version={version}/dataset.rds",
  name = "Chronos_phase1_proteomics_somalogic_analysis_ready_three_tibble",
  version = 29650
) |>
  readr::read_rds()
somalogic_table <- somalogic_tibble$meta_data |>
  dplyr::inner_join(
    somalogic_tibble$prot_data,
    by = "ID"
  )
dim(somalogic_table)

somalogic_results <- glue::glue(
  "{get_config()$input_dir}/name={name}/version={version}/dataset.parquet",
  name = "chronos_somalogic_merged-pd_master_summary_2",
  version = 29828
) |>
  arrow::open_dataset() |>
  dplyr::collect()

somalogic_JM <- glue::glue(
  "{get_config()$input_dir}/name={name}/version={version}/dataset.parquet",
  name = "chronos_somalogic_joint_model-pd_master_summary",
  version = 29819
) |>
  arrow::open_dataset() |>
  dplyr::collect()

#####
# olink ####
#####
olink_tibble <- glue::glue(
  "{get_config()$input_dir}/name={name}/version={version}/dataset.rds",
  name = "Chronos_phase1_proteomics_olink_analysis_ready_three_tibble",
  version = 29465
) |>
  readr::read_rds()

olink_table <- olink_tibble$meta_data |>
  dplyr::inner_join(
    olink_tibble$prot_data,
    by = "ID"
  )

olink_results <- glue::glue(
  "{get_config()$input_dir}/name={name}/version={version}/dataset.parquet",
  name = "chronos_olink_merged-pd_master_summary_2",
  version = 29829
) |>
  arrow::open_dataset() |>
  dplyr::collect()

olink_JM <- glue::glue(
  "{get_config()$input_dir}/name={name}/version={version}/dataset.parquet",
  name = "chronos_olink_joint_model-pd_master_summary",
  version = 29820
) |>
  arrow::open_dataset() |>
  dplyr::collect()

#####
# alamar ####
#####
alamar_tibble <- glue::glue(
  "{get_config()$input_dir}/name={name}/version={version}/dataset.rds",
  name = "Chronos_phase1_proteomics_alamar_analysis_ready_three_tibble",
  version = 29649
) |>
  readr::read_rds()

alamar_table <- alamar_tibble$meta_data |>
  dplyr::inner_join(
    alamar_tibble$prot_data,
    by = "ID"
  )
alamar_results <- glue::glue(
  "{get_config()$input_dir}/name={name}/version={version}/dataset.parquet",
  name = "chronos_alamar_merged-pd_master_summary_2",
  version = 29830
) |>
  arrow::open_dataset() |>
  dplyr::collect()
alamar_JM <- glue::glue(
  "{get_config()$input_dir}/name={name}/version={version}/dataset.parquet",
  name = "chronos_alamar_joint_model-pd_master_summary-pd_master_summary_2",
  version = 29821
) |>
  arrow::open_dataset() |>
  dplyr::collect()


#####
# chronic condition ICD10 ####
#####

infile_ccir_icd10_3digits <- glue::glue("{get_config()$input_dir}/name=ccir_icd10_3digits/version=none/dataset.parquet")

phase1_chronic <- infile_ccir_icd10_3digits |>
  arrow::open_dataset() |>
  dplyr::collect() |>
  dplyr::select(
    PATIENT_ID = patient_id,
    DIAGNOSIS_DATE = min_diagnosis_date,
    diagnosis_code_3 = diagnosis_code_3digits
  )

# OR on chronic 3-digit ICD10 ####
df_codes <- phase1_critical_variables |>
  dplyr::distinct(patient_id, group) |>
  dplyr::left_join(
    phase1_chronic,
    by = c("patient_id" = "PATIENT_ID")
  )

totals <- df_codes |>
  dplyr::distinct(patient_id, group) |>
  dplyr::count(group) |>
  tidyr::pivot_wider(names_from = group, values_from = n)

n_pd <- totals$pd # 348
n_control <- totals$control # 348

counts <- df_codes |>
  dplyr::count(diagnosis_code_3, group) |>
  tidyr::pivot_wider(names_from = group, values_from = n, values_fill = 0L) |>
  dplyr::rename(
    a = pd, # pd with code
    b = control
  ) # control with code
head(counts)

results <- counts |>
  dplyr::mutate(
    # patients without the code
    c = n_pd - a,
    d = n_control - b,
    # run Fisher’s exact test:
    fisher = purrr::map2(a, b, ~ fisher.test(matrix(c(.x, .y, n_pd - .x, n_control - .y),
      nrow = 2,
      byrow = TRUE
    ))),
    # extract OR, CIs, p
    OR = purrr::map_dbl(fisher, ~ .$estimate),
    CI_lower = purrr::map_dbl(fisher, ~ .$conf.int[1]),
    CI_upper = purrr::map_dbl(fisher, ~ .$conf.int[2]),
    p_value = purrr::map_dbl(fisher, ~ .$p.value)
  ) |>
  dplyr::select(diagnosis_code_3, OR, CI_lower, CI_upper, p_value, a, b) |>
  dplyr::rename(
    n_pd = a,
    n_control = b
  ) |>
  dplyr::mutate(n_total = n_pd + n_control) |>
  dplyr::arrange(p_value)


OR_phase1_chronic <- results |>
  dplyr::mutate(
    adjp = ifelse(
      diagnosis_code_3 %in% c("G20", "G21"), # adjust without having G20, and G21
      NA_real_,
      p.adjust(p_value[diagnosis_code_3 %nin% c("G20", "G21")], method = "BH")[match(diagnosis_code_3, diagnosis_code_3[diagnosis_code_3 %nin% c("G20", "G21")])]
    )
  ) |>
  dplyr::filter(!is.na(diagnosis_code_3))
dim(OR_phase1_chronic)

save_versioned_dataset(OR_phase1_chronic, "OR_phase1_chronic_updated0120")

OR_phase1_chronic_ref <- infile_ccir_icd10_3digits |>
  arrow::open_dataset() |>
  dplyr::collect() |>
  dplyr::filter(diagnosis_code_3digits %in% OR_phase1_chronic$diagnosis_code_3) |>
  dplyr::distinct(diagnosis_code_3digits, icd10_description_3digits) |>
  dplyr::mutate(DIAGNOSIS_CODE = diagnosis_code_3digits) |>
  dplyr::rename(ICD10_description = icd10_description_3digits,
                diagnosis_code_3 = diagnosis_code_3digits) |>
  dplyr::select(DIAGNOSIS_CODE, ICD10_description, diagnosis_code_3)

save_versioned_dataset(OR_phase1_chronic_ref, "OR_phase1_chronic_ref")

significant_chronic <- OR_phase1_chronic |>
  dplyr::filter(adjp < 0.05)

significant_chronic <- significant_chronic |>
  dplyr::inner_join(
    OR_phase1_chronic_ref,
    by = "diagnosis_code_3"
  )
dim(significant_chronic)

# get the first comorbidity ####
first_comorbidity <- phase1_chronic |>
  dplyr::filter(diagnosis_code_3 %in% c(significant_chronic$diagnosis_code_3, "G20")) |> # filtering on chronic conditions only
  dplyr::group_by(PATIENT_ID, diagnosis_code_3) |>
  dplyr::filter(DIAGNOSIS_DATE == min(DIAGNOSIS_DATE)) |> # get the first date
  dplyr::ungroup() |>
  dplyr::distinct(PATIENT_ID, DIAGNOSIS_DATE, diagnosis_code_3)
dim(first_comorbidity)
length(unique(first_comorbidity$PATIENT_ID))
length(unique(first_comorbidity$diagnosis_code_3))
# get the first non-PD related conditions ####
# Focus on non-PD related conditions
non_significant_chronic <- OR_phase1_chronic |>
  dplyr::filter(adjp >= 0.05) |>
  dplyr::filter(n_total >= 30)
dim(non_significant_chronic)
first_condition <- phase1_chronic |>
  dplyr::filter(diagnosis_code_3 %in% c(non_significant_chronic$diagnosis_code_3)) |> # filtering on chronic conditions only
  dplyr::group_by(PATIENT_ID, diagnosis_code_3) |>
  dplyr::filter(DIAGNOSIS_DATE == min(DIAGNOSIS_DATE)) |> # get the first date
  dplyr::ungroup() |>
  dplyr::distinct(PATIENT_ID, DIAGNOSIS_DATE, diagnosis_code_3)
dim(first_condition)
length(unique(first_condition$PATIENT_ID))
length(unique(first_condition$diagnosis_code_3))

# prepare survival table ####
#####
# read phase 1 cohort ICD10 converted claim records:
patient_claims_converted <- glue::glue("{get_config()$input_dir}/name=patient_claims_converted/version=none/dataset.parquet")

claim_summary_phase1 <- patient_claims_converted |>
  arrow::open_dataset() |>
  dplyr::collect() |>
  dplyr::group_by(patient_id) |>
  dplyr::summarise(
    first_date_id = min(diagnosis_date),
    last_date_id = max(diagnosis_date)
  ) |>
  dplyr::rename(PATIENT_ID = patient_id)

dim(claim_summary_phase1)

first_comorbidity_summary <- claim_summary_phase1 |>
  dplyr::left_join(
    first_comorbidity,
    by = "PATIENT_ID"
  )
dim(first_comorbidity_summary)
all_codes <- unique(first_comorbidity$diagnosis_code_3) # 51, inlcude 20
first_comorbidity_summary_completed <- claim_summary_phase1 |>
  tidyr::crossing(diagnosis_code_3 = all_codes) |>
  dplyr::left_join(
    first_comorbidity_summary,
    by = c("PATIENT_ID", "first_date_id", "last_date_id", "diagnosis_code_3")
  ) |>
  dplyr::arrange(PATIENT_ID, diagnosis_code_3, DIAGNOSIS_DATE) |>
  dplyr::left_join(
    phase1_critical_variables |>
      dplyr::distinct(patient_id, group, birth_year),
    by = c("PATIENT_ID" = "patient_id")
  ) |>
  dplyr::group_by(PATIENT_ID, diagnosis_code_3) |>
  dplyr::mutate(
    Event_date = as.Date(ifelse(is.na(DIAGNOSIS_DATE), last_date_id, DIAGNOSIS_DATE)), # one patient with NA
    birth_year = as.Date(paste0(birth_year, "-01-01"))
  ) |>
  dplyr::ungroup()
dim(first_comorbidity_summary_completed)
first_condition_summary <- claim_summary_phase1 |>
  dplyr::left_join(
    first_condition,
    by = "PATIENT_ID"
  )
dim(first_condition_summary)
all_codes <- unique(first_condition$diagnosis_code_3) # 52
first_condition_summary_completed <- claim_summary_phase1 |>
  tidyr::crossing(diagnosis_code_3 = all_codes) |>
  dplyr::left_join(
    first_condition_summary,
    by = c("PATIENT_ID", "first_date_id", "last_date_id", "diagnosis_code_3")
  ) |>
  dplyr::arrange(PATIENT_ID, diagnosis_code_3, DIAGNOSIS_DATE) |>
  dplyr::left_join(
    phase1_critical_variables |>
      dplyr::distinct(patient_id, group, birth_year),
    by = c("PATIENT_ID" = "patient_id")
  ) |>
  dplyr::group_by(PATIENT_ID, diagnosis_code_3) |>
  dplyr::mutate(
    Event_date = as.Date(ifelse(is.na(DIAGNOSIS_DATE), last_date_id, DIAGNOSIS_DATE)), # one patient with NA
    birth_year = as.Date(paste0(birth_year, "-01-01"))
  ) |>
  dplyr::ungroup()
dim(first_condition_summary_completed)
#####
# merge with sample records ####
sample_comorbidity_combined <- first_comorbidity_summary_completed |>
  dplyr::select(PATIENT_ID, group, diagnosis_code_3, DIAGNOSIS_DATE, Event_date, first_date_id, birth_year) |>
  dplyr::left_join(
    phase1_critical_variables |>
      dplyr::distinct(patient_id, gender, race, sample_id, sample_collection_date, recent_plasma_donations, cluster),
    by = c("PATIENT_ID" = "patient_id")
  ) |>
  dplyr::filter(sample_collection_date < Event_date) |> # only get samples before event
  dplyr::mutate(
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
    )) + 1) / 365.25
  )
dim(sample_comorbidity_combined)
sample_condition_combined <- first_condition_summary_completed |>
  dplyr::select(PATIENT_ID, group, diagnosis_code_3, DIAGNOSIS_DATE, Event_date, first_date_id, birth_year) |>
  dplyr::left_join(
    phase1_critical_variables |>
      dplyr::distinct(patient_id, gender, race, sample_id, sample_collection_date, recent_plasma_donations, cluster),
    by = c("PATIENT_ID" = "patient_id")
  ) |>
  dplyr::filter(sample_collection_date < Event_date) |> # only get samples before event
  dplyr::mutate(
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
    )) + 1) / 365.25
  )
dim(sample_condition_combined)


# list of biomarkers to test ####
risk_combined <- somalogic_JM |>
  dplyr::filter(rhat < 1.05) |>
  dplyr::select(probe_id, uniprot_id, gene_symbol, estimate, p_val, q_val, risk, risk_lower, risk_upper) |>
  dplyr::mutate(platform = "SomaScanV5") |>
  dplyr::bind_rows(olink_JM |>
    dplyr::filter(rhat < 1.05) |>
    dplyr::select(probe_id, uniprot_id, gene_symbol, estimate, p_val, q_val, risk, risk_lower, risk_upper) |>
    dplyr::mutate(platform = "Olink_HT")) |>
  dplyr::bind_rows(alamar_JM |>
    dplyr::filter(rhat < 1.05) |>
    dplyr::select(probe_id, uniprot_id, gene_symbol, estimate, p_val, q_val, risk, risk_lower, risk_upper) |>
    dplyr::mutate(platform = "Alamarbio_CNS")) |>
  dplyr::group_by(uniprot_id) |>
  dplyr::filter(p_val == min(p_val)) |>
  dplyr::ungroup() |>
  dplyr::mutate(gene_symbol = reorder(gene_symbol, risk)) |>
  dplyr::filter(p_val < 0.005)
dim(risk_combined)
early_combined <- somalogic_results |>
  dplyr::filter(q_val_pre < 0.05) |>
  dplyr::group_by(uniprot_id) |>
  dplyr::filter(p_val_pre == min(p_val_pre)) |>
  dplyr::ungroup() |>
  dplyr::select(
    probe_id, uniprot_id, gene_symbol,
    log2_fc_pre, p_val_pre, cohens_d_pre, q_val_pre
  ) |>
  dplyr::mutate(platform = "SomaScanV5") |>
  dplyr::bind_rows(
    olink_results |>
      dplyr::filter(q_val_pre < 0.05) |>
      dplyr::group_by(uniprot_id) |>
      dplyr::filter(p_val_pre == min(p_val_pre)) |>
      dplyr::ungroup() |>
      dplyr::select(
        probe_id, uniprot_id, gene_symbol,
        log2_fc_pre, p_val_pre, cohens_d_pre, q_val_pre
      ) |>
      dplyr::mutate(platform = "Olink_HT")
  ) |>
  dplyr::bind_rows(
    alamar_results |>
      dplyr::filter(q_val_pre < 0.05) |>
      dplyr::group_by(uniprot_id) |>
      dplyr::filter(p_val_pre == min(p_val_pre)) |>
      dplyr::ungroup() |>
      dplyr::select(
        probe_id, uniprot_id, gene_symbol,
        log2_fc_pre, p_val_pre, cohens_d_pre, q_val_pre
      ) |>
      dplyr::mutate(platform = "Alamarbio_CNS")
  )
dim(early_combined)
protein_list <- risk_combined |>
  dplyr::distinct(probe_id, uniprot_id, gene_symbol, platform) |>
  dplyr::bind_rows(
    early_combined |>
      dplyr::distinct(probe_id, uniprot_id, gene_symbol, platform)
  ) |>
  dplyr::distinct()
dim(protein_list)
save_versioned_dataset(protein_list, "protein_list")


# prepare JM table ####

#####
# somalogic:
#####
JM_lme_somalogic <- sample_comorbidity_combined |>
  dplyr::inner_join(
    somalogic_table |>
      dplyr::select(sample_id, tidyselect::all_of(protein_list$probe_id[which(protein_list$platform == "SomaScanV5")])),
    by = "sample_id"
  )
dim(JM_lme_somalogic)
JM_surv_somalogic <- JM_lme_somalogic |>
  dplyr::distinct(
    PATIENT_ID, diagnosis_code_3, age_event, age_entry, Event_date, DIAGNOSIS_DATE,
    group, gender, race, birth_year
  ) |>
  dplyr::mutate(
    Time_age = age_event,
    Status = ifelse(is.na(DIAGNOSIS_DATE), 0, 1)
  )
dim(JM_surv_somalogic)
save_versioned_dataset(JM_lme_somalogic, "JM_lme_somalogic")
save_versioned_dataset(JM_surv_somalogic, "JM_surv_somalogic")

#####
# olink:
#####
JM_lme_olink <- sample_comorbidity_combined |>
  dplyr::inner_join(
    olink_table |>
      dplyr::select(sample_id, tidyselect::all_of(protein_list$probe_id[which(protein_list$platform == "Olink_HT")])),
    by = "sample_id"
  )
dim(JM_lme_olink)
JM_surv_olink <- JM_lme_olink |>
  dplyr::distinct(
    PATIENT_ID, diagnosis_code_3, age_event, age_entry, Event_date, DIAGNOSIS_DATE,
    group, gender, race, birth_year
  ) |>
  dplyr::mutate(
    Time_age = age_event,
    Status = ifelse(is.na(DIAGNOSIS_DATE), 0, 1)
  )
dim(JM_surv_olink)
save_versioned_dataset(JM_lme_olink, "JM_lme_olink")
save_versioned_dataset(JM_surv_olink, "JM_surv_olink")

#####
# alamar:
#####
JM_lme_alamar <- sample_comorbidity_combined |>
  dplyr::inner_join(
    alamar_table |>
      dplyr::select(sample_id, tidyselect::all_of(protein_list$probe_id[which(protein_list$platform == "Alamarbio_CNS")])),
    by = "sample_id"
  )
dim(JM_lme_alamar)
JM_surv_alamar <- JM_lme_alamar |>
  dplyr::distinct(
    PATIENT_ID, diagnosis_code_3, age_event, age_entry, Event_date, DIAGNOSIS_DATE,
    group, gender, race, birth_year
  ) |>
  dplyr::mutate(
    Time_age = age_event,
    Status = ifelse(is.na(DIAGNOSIS_DATE), 0, 1)
  )
dim(JM_surv_alamar)
save_versioned_dataset(JM_lme_alamar, "JM_lme_alamar")
save_versioned_dataset(JM_surv_alamar, "JM_surv_alamar")

#####
# somalogic:
#####
JMcondition_lme_somalogic <- sample_condition_combined |>
  dplyr::inner_join(
    somalogic_table |>
      dplyr::select(sample_id, tidyselect::all_of(protein_list$probe_id[which(protein_list$platform == "SomaScanV5")])),
    by = "sample_id"
  )
dim(JMcondition_lme_somalogic)
JMcondition_surv_somalogic <- JMcondition_lme_somalogic |>
  dplyr::distinct(
    PATIENT_ID, diagnosis_code_3, age_event, age_entry, Event_date, DIAGNOSIS_DATE,
    group, gender, race, birth_year
  ) |>
  dplyr::mutate(
    Time_age = age_event,
    Status = ifelse(is.na(DIAGNOSIS_DATE), 0, 1)
  )
dim(JMcondition_surv_somalogic)
save_versioned_dataset(JMcondition_lme_somalogic, "JMcondition_lme_somalogic")
save_versioned_dataset(JMcondition_surv_somalogic, "JMcondition_surv_somalogic")

#####
# olink:
#####
JMcondition_lme_olink <- sample_condition_combined |>
  dplyr::inner_join(
    olink_table |>
      dplyr::select(sample_id, tidyselect::all_of(protein_list$probe_id[which(protein_list$platform == "Olink_HT")])),
    by = "sample_id"
  )
dim(JMcondition_lme_olink)
JMcondition_surv_olink <- JMcondition_lme_olink |>
  dplyr::distinct(
    PATIENT_ID, diagnosis_code_3, age_event, age_entry, Event_date, DIAGNOSIS_DATE,
    group, gender, race, birth_year
  ) |>
  dplyr::mutate(
    Time_age = age_event,
    Status = ifelse(is.na(DIAGNOSIS_DATE), 0, 1)
  )
dim(JMcondition_surv_olink)
save_versioned_dataset(JMcondition_lme_olink, "JMcondition_lme_olink")
save_versioned_dataset(JMcondition_surv_olink, "JMcondition_surv_olink")

#####
# alamar:
#####
JMcondition_lme_alamar <- sample_condition_combined |>
  dplyr::inner_join(
    alamar_table |>
      dplyr::select(sample_id, tidyselect::all_of(protein_list$probe_id[which(protein_list$platform == "Alamarbio_CNS")])),
    by = "sample_id"
  )
dim(JMcondition_lme_alamar)
JMcondition_surv_alamar <- JMcondition_lme_alamar |>
  dplyr::distinct(
    PATIENT_ID, diagnosis_code_3, age_event, age_entry, Event_date, DIAGNOSIS_DATE,
    group, gender, race, birth_year
  ) |>
  dplyr::mutate(
    Time_age = age_event,
    Status = ifelse(is.na(DIAGNOSIS_DATE), 0, 1)
  )
dim(JMcondition_surv_alamar)
save_versioned_dataset(JMcondition_lme_alamar, "JMcondition_lme_alamar")
save_versioned_dataset(JMcondition_surv_alamar, "JMcondition_surv_alamar")

# Summary ####
cli::cli_h2("Data preparation complete!")
cli::cli_alert_info(glue::glue("All output files saved to: {data_prep_base}/"))
cli::cli_alert_info("Files are organized in versioned format: name=<dataset>/version=<version>/dataset.parquet")
cli::cli_alert_info("Note: These files should be published to the data publishing repository for use by modeling scripts")
