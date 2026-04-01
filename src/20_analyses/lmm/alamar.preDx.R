## SPDX-FileCopyrightText: 2026 Alkahest, Inc. <chronospd_data_access@alkahest.com>
##
## SPDX-License-Identifier: MIT

here::i_am("src/20_analyses/lmm/alamar.preDx.R")

# screen model ####
input_alamar_dataset_file <- glue::glue("{get_config()$input_dir}/name=Chronos_phase1_proteomics_alamar_analysis_ready_three_tibble/version=29464/dataset.rds")

output_intermediate_result_file <- glue::glue("{get_config()$output_dir}/lmm/intermediate/alamar.preDx.screening.rds")
output_result_file <- glue::glue("{get_config()$output_dir}/lmm/alamar.preDx.screening.rds")

cores <- parallelly::availableCores()

dir.create(dirname(output_intermediate_result_file), showWarnings = FALSE, recursive = TRUE)
dir.create(dirname(output_result_file), showWarnings = FALSE, recursive = TRUE)

alamar_tibble <- readRDS(input_alamar_dataset_file)

# preDx model ####
ID_pre <- alamar_tibble$meta_data |>
  dplyr::filter(years_since_pd_onset < 0)

alamar_tibble_pre <- list(
  prot_data = alamar_tibble$prot_data |>
    dplyr::filter(ID %in% ID_pre$ID) |>
    tibble::as_tibble(),
  meta_data = ID_pre |>
    tibble::as_tibble(),
  prot_annot = alamar_tibble$prot_annot |>
    tibble::as_tibble()
)

fixed_effects_list <- c(
  "years_since_pd_onset*group + gender + age_at_baseline + recent_plasma_donations + cluster",
  "years_since_pd_onset*group + gender + age_at_baseline + recent_plasma_donations + cluster + race"
)

random_effects_list <- c(
  "1 | patient_id",
  "1 | patient_id + 1 | center",
  "1 + years_since_pd_onset | patient_id",
  "1 | patient_id + (0 + years_since_pd_onset | patient_id)"
)

res <- pbmcapply::pbmclapply(alamar_tibble_pre$prot_annot |>
                               dplyr::pull(ID), # Iterate through protein IDs
                             compare_multiple_models_low_mem_wrapper, ### <- The function being ran
                             three_tibble_list = alamar_tibble_pre,
                             fixed_effects_list = fixed_effects_list,
                             random_effects_list = random_effects_list,
                             mc.cores = cores)

saveRDS(res, output_intermediate_result_file)

grouppd_table <- res |>
  purrr::map(function(x) {

    ttab <- x$best_model_summary$tTable
    prot <- x$best_model_summary$protein

    if (base::is.null(ttab) || base::is.null(prot)) return(NULL)

    id_val <- prot$ID[1]

    tt <- tibble::as_tibble(ttab) |>
      dplyr::filter(coef.names == "grouppd")

    if (base::nrow(tt) == 0) return(NULL)

    wanted_cols <- c("coef.names","Estimate","std.err","DF","t.value","p.value","cohensD")

    tt <- tt[, base::intersect(base::names(tt), wanted_cols), drop = FALSE]

    tt$ID <- id_val
    tt
  }) |>
  purrr::compact() |>
  dplyr::bind_rows() |>
  dplyr::mutate(adj.p = stats::p.adjust(p.value, method = "BH")) |>
  dplyr::arrange(p.value) |>
  dplyr::relocate(ID)

grouppd_table |>
  dplyr::filter(adj.p < 0.05) # 7

saveRDS(grouppd_table, output_result_file)
