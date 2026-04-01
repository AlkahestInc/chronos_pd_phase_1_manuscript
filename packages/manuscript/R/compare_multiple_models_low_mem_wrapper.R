## SPDX-FileCopyrightText: 2026 Alkahest, Inc. <chronospd_data_access@alkahest.com>
##
## SPDX-License-Identifier: MIT

#' compare_multiple_models_low_mem_wrapper
#'
#' Low-memory wrapper around `compare_multiple_models`. Re-fits the best model with
#' `lmerTest` and `REML = TRUE` for inference and adds diagnostics & R2.
#'
#' @param protein_id Character.
#' @param three_tibble_list See `derive_model_data()`.
#' @param ... Arguments forwarded to `compare_multiple_models`.
#' @return A list with model_summaries and best_model_summary.
#' @export
compare_multiple_models_low_mem_wrapper <- function(protein_id, three_tibble_list, ...) {
  model_results <- compare_multiple_models(protein_id, three_tibble_list, package = "lmerTest", REML = FALSE, ...)
  if (is.null(model_results$best_model)) {
    return(list(model_summaries = model_results$models, best_model_summary = NULL))
  }

  # Robust refit with REML=TRUE (no stats::update)
  fit <- explain_by_LME(
    protein_id,
    three_tibble_list,
    fixed_effects = model_results$best_model_fixed,
    random_effects = model_results$best_model_random,
    package = "lmerTest",
    REML = TRUE
  )

  protein_id_chr <- as.character(stats::formula(fit))[2]
  prot_annot <- three_tibble_list$prot_annot |> dplyr::filter(.data$ID == protein_id_chr)

  tTable <- extract_model_summary_statistics(fit) |>
    dplyr::rename(std.err = "Std.Error")

  # Diagnostics + R2
  pTable <- extract_model_diagnostic_statistics(fit)
  r2 <- tryCatch(performance::r2(fit), error = function(e) NULL)
  if (!is.null(r2)) {
    if (isTRUE(all(c("R2_conditional", "R2_marginal") %in% names(r2)))) {
      pTable$R2_conditional <- r2[["R2_conditional"]]
      pTable$R2_marginal <- r2[["R2_marginal"]]
    } else if (is.data.frame(r2) && all(c("R2_conditional","R2_marginal") %in% colnames(r2))) {
      pTable$R2_conditional <- r2$R2_conditional
      pTable$R2_marginal <- r2$R2_marginal
    } else {
      pTable$R2_conditional <- NA_real_
      pTable$R2_marginal <- NA_real_
    }
  } else {
    pTable$R2_conditional <- NA_real_
    pTable$R2_marginal <- NA_real_
  }

  list(
    model_summaries = model_results$models,
    best_model_summary = list(
      formula = paste(deparse(stats::formula(model_results$best_model)), collapse = ""),
      BIC = model_results$best_bic,
      protein = prot_annot,
      tTable = tTable,
      pTable = pTable
    )
  )
}
