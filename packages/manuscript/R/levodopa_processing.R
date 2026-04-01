## SPDX-FileCopyrightText: 2026 Alkahest, Inc. <chronospd_data_access@alkahest.com>
##
## SPDX-License-Identifier: MIT

#' Align post-PD samples to medications and calculating the levodopa effect
#'
#' @param dataprot Raw protein tables (rows: samples; columns: probes)
#' @param dataother Sample annotation
#' @param prot_annot Probe annotation
#' @param patient_meds_df RWD medication records
#' @param phase1_donors_df Donor annotation
#' @param pd_meds_fun levodopa containing NDCs
#'
#' @return A ggplot2 object

levodopa_processing <- function(
    dataprot,
    dataother,
    prot_annot,
    patient_meds_df,
    phase1_donors_df,
    pd_meds_fun,
    time_of_medication_effect = -30,
    time_of_exclusion_medication_effect = -365,
    other_ongoing_medication_window = 90
) {
  # medication-sample post-PD alignment
  patient_meds_selected <- patient_meds_df |>
    dplyr::right_join(
      dataother |>
        dplyr::filter(
          years_since_pd_onset >= 0,
          group == "pd"
        ) |>
        dplyr::select(patient_id, aliquot_barcode, sample_collection_date),
      by = "patient_id",
      relationship = "many-to-many"
    ) |>
    dplyr::mutate(time_from_sample = as.numeric(date_of_service - sample_collection_date)) |>
    dplyr::left_join(pd_meds_fun, by = "national_drug_code") |>
    dplyr::mutate(MED_GROUP = dplyr::if_else(is.na(MED_GROUP), "other_medication", MED_GROUP))
  # if the time_from_sample is negative that means that the medication happened before the sample

  # 1 selecting the levodopa+ samples (samples with prior levodopa)
  with_levodopa_samples <- patient_meds_selected |>
    dplyr::filter(MED_GROUP != "other_medication") |>
    dplyr::filter(
      time_from_sample > time_of_medication_effect,
      time_from_sample < 0
    ) |>
    dplyr::mutate(group_name = "with_levodopa") |>
    dplyr::select(aliquot_barcode, group_name) |>
    dplyr::distinct()

  # 2 selecting the levodopa- samples (samples with prior levodopa)
  # 2A selecting samples where we see prior levodopa in a wide temporal window
  levodopa_suspected_samples <- patient_meds_selected |>
    dplyr::filter(MED_GROUP != "other_medication") |>
    dplyr::filter(
      time_from_sample > time_of_exclusion_medication_effect,
      time_from_sample < 0
    ) |>
    dplyr::select(aliquot_barcode) |>
    dplyr::distinct()

  # 2B selecting those samples where we don't see levodopa in a wide window (from 2A) and has OTHER medications before and after the samples.
  without_levodopa_samples <- patient_meds_selected |>
    dplyr::filter(
      !(aliquot_barcode %in% levodopa_suspected_samples$aliquot_barcode),
    ) |>
    dplyr::group_by(aliquot_barcode) |>
    dplyr::summarise(
      has_neg_meds = any(time_from_sample >= (-1 * other_ongoing_medication_window) & time_from_sample <= 0),
      has_pos_meds = any(time_from_sample >= 0 & time_from_sample <= other_ongoing_medication_window)
    ) |>
    dplyr::filter(has_neg_meds & has_pos_meds) |>
    dplyr::mutate(group_name = "without_levodopa") |>
    dplyr::select(aliquot_barcode, group_name) |>
    dplyr::distinct()

  rwd_case_samples <- dplyr::bind_rows(
    with_levodopa_samples |> dplyr::select(aliquot_barcode, group_name),
    without_levodopa_samples |> dplyr::select(aliquot_barcode, group_name)
  )

  # DATA PREPARATION for LMM - Levodopa
  subset_meta <- dataother |>
    dplyr::inner_join(rwd_case_samples, by = "aliquot_barcode")

  ### donors with "without_levodopa" and "with_levodopa" samples were assigned to the group where the majority of their samples came from (only 4 donors)
  test_multiple <- subset_meta |>
    dplyr::group_by(patient_id, group_name) |>
    dplyr::summarise(num_group = dplyr::n(), .groups = "drop")

  patients_both <- test_multiple |>
    dplyr::group_by(patient_id) |>
    dplyr::filter(dplyr::n() == 2) |>
    dplyr::ungroup()

  patients_remove <- patients_both |>
    dplyr::group_by(patient_id) |>
    dplyr::slice_min(num_group, n = 1, with_ties = FALSE) |>
    dplyr::ungroup()

  meta <- subset_meta |>
    dplyr::anti_join(
      patients_remove |> dplyr::select(patient_id, group_name),
      by = c("patient_id", "group_name")
    )
  ###

  prot_data <- dataprot |>
    dplyr::filter(aliquot_barcode %in% meta$aliquot_barcode) |>
    dplyr::select(aliquot_barcode, dplyr::all_of(prot_annot$protein_id))

  protein_ids <- prot_annot$protein_id

  ## LMM
  cores <- parallelly::availableCores()
  future::plan(future::multisession, workers = cores)

  compiled_tt <- furrr::future_map(protein_ids, function(pid) {
    df <- meta |>
      dplyr::left_join(
        prot_data |> dplyr::select(aliquot_barcode, !!pid),
        by = "aliquot_barcode"
      ) |>
      dplyr::rename(y = !!pid) |>
      dplyr::filter(!is.na(y))
    df$group_name <- factor(df$group_name, levels = c("without_levodopa", "with_levodopa"))
    # LMM fixed + random intercept
    f <- stats::as.formula(
      "y ~ years_since_pd_onset + group_name +
            gender + age_at_baseline +
            recent_plasma_donations + cluster +
            (1 | patient_id)"
    )

    fit <- lmerTest::lmer(f,
      data = df,
      REML = TRUE,
      control = lme4::lmerControl(
        optimizer = "bobyqa",
        optCtrl   = list(maxfun = 30000)
      )
    )

    # Extract fixed effects
    coefmat <- summary(fit)$coefficients
    tTable <- tibble::as_tibble(coefficients(summary(fit)), rownames = "coef.names") |>
      dplyr::rename(t.value = "t value", p.value = "Pr(>|t|)", std.err = "Std. Error")

    tt <- data.frame(
      term = rownames(coefmat),
      estimate = coefmat[, "Estimate"],
      std_error = coefmat[, "Std. Error"],
      df = coefmat[, "df"],
      t_value = coefmat[, "t value"],
      p_value = coefmat[, "Pr(>|t|)"],
      stringsAsFactors = FALSE
    )
    rownames(tt) <- NULL

    # Convergence
    opt <- fit@optinfo$conv$opt
    msgs <- fit@optinfo$conv$lme4$messages
    opt <- fit@optinfo$conv$opt

    # convergence messages (may include harmless boundary warnings)
    msgs <- fit@optinfo$conv$lme4$messages

    # 1. TRUE convergence means: optimizer succeeded (opt == 0)
    tt$converged <- isTRUE(opt == 0)

    # 2. Singularity is independent of convergence:
    tt$singular <- lme4::isSingular(fit, tol = 1e-4)

    cohensd <- effectsize::t_to_d(tTable[["t.value"]], tTable[["df"]])
    tt$cohensD <- cohensd$d

    tt$protein_id <- pid

    tt
  }) |>
    dplyr::bind_rows()

  postonset_rwdplus_table <- compiled_tt |>
    dplyr::filter(grepl("^group_name", term)) |>
    dplyr::mutate(p.value = p_value) |>
    dplyr::mutate(adj = p.adjust(p.value, method = "BH")) |>
    dplyr::arrange(adj) |>
    dplyr::mutate(
      neg_log10_p = -log10(p.value)
    )

  # levodopa LMM results
  lmm_results <- postonset_rwdplus_table |>
    dplyr::select(protein_id, cohensD, p.value, neg_log10_p, estimate, std_error, t_value, df) |>
    dplyr::inner_join(
      prot_annot |>
        dplyr::select(protein_id, gene_symbol),
      by = "protein_id"
    ) |>
    dplyr::relocate(gene_symbol, .after = protein_id)

  lmm_results
}
