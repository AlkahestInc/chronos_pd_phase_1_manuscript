## SPDX-FileCopyrightText: 2026 Alkahest, Inc. <chronospd_data_access@alkahest.com>
##
## SPDX-License-Identifier: MIT

#' Analyze a single protein ratio using logistic regression
#'
#' Fits a logistic regression model for a single protein ratio
#' (log2 intensity difference) against PD status, adjusting for
#' specified covariates. Returns a one-row tibble with model results
#' or an empty-result template on failure.
#'
#' @param analyte_1 Character. Column name for protein 1 intensity.
#' @param analyte_2 Character. Column name for protein 2 intensity.
#' @param protein_pair Character. Identifier for the protein pair.
#' @param prot_1_uniprot Character. UniProt ID for protein 1.
#' @param prot_2_uniprot Character. UniProt ID for protein 2.
#' @param prot_1_symbol Character. Gene symbol for protein 1.
#' @param prot_2_symbol Character. Gene symbol for protein 2.
#' @param protein_ratio_metadata Data frame with sample metadata.
#' @param protein_intensity_log2 Data frame of log2 protein intensities.
#' @param other_predictors Character vector of covariate column names. Default to c("age_at_visit", "sex")
#'
#' @return One-row tibble with regression results or empty-result template.
analyze_single_ratio <- function(
  analyte_1,
  analyte_2,
  protein_pair,
  prot_1_uniprot,
  prot_2_uniprot,
  prot_1_symbol,
  prot_2_symbol,
  protein_ratio_metadata,
  protein_intensity_log2,
  other_predictors = c("age_at_visit", "sex")
) {
  # Input validation
  if (!(is.character(protein_pair) && length(protein_pair) == 1)) {
    stop(paste0("Parameter protein_pair must be a single string: ", protein_pair))
  }

  # Empty result template
  empty_result <- tibble::tibble(
    protein_pair = protein_pair,
    mapped_analyte_1 = ifelse(length(analyte_1) == 0, NA_character_, analyte_1),
    mapped_analyte_2 = ifelse(length(analyte_2) == 0, NA_character_, analyte_2),
    prot_1_uniprot = ifelse(length(prot_1_uniprot) == 0, NA_character_, prot_1_uniprot),
    prot_2_uniprot = ifelse(length(prot_2_uniprot) == 0, NA_character_, prot_2_uniprot),
    gene_symbol1 = ifelse(length(prot_1_symbol) == 0, NA_character_, prot_1_symbol),
    gene_symbol2 = ifelse(length(prot_2_symbol) == 0, NA_character_, prot_2_symbol),
    estimate = NA_real_,
    or_lower_ci = NA_real_,
    or_upper_ci = NA_real_,
    p_value = NA_real_,
    n_obs = NA_integer_,
    converged = FALSE
  )

  # Check analyte existence
  if (!(is.character(analyte_1) && length(analyte_1) == 1)) {
    # stop(paste0("Parameter analyte_1 must be a single string: ", analyte_1))
    return(empty_result)
  }
  if (!(is.character(analyte_2) && length(analyte_2) == 1)) {
    # stop(paste0("Parameter analyte_2 must be a single string: ", analyte_2))
    return(empty_result)
  }
  if (!(analyte_1 %in% annotations$column_name && analyte_2 %in% annotations$column_name)) {
    return(empty_result)
  }

  result <- tryCatch(
    {
      protein_intensity_log2_by_prot_1 <- protein_intensity_log2 |>
        select(prot_1_symbol = all_of(analyte_1))
      protein_intensity_log2_by_prot_2 <- protein_intensity_log2 |>
        select(prot_2_symbol = all_of(analyte_2))

      ratio_name <- paste0(prot_1_symbol, "_vs_", prot_2_symbol)
      ratio_df <- protein_intensity_log2 |>
        select(sample_id, prot_1_symbol = all_of(analyte_1), prot_2_symbol = all_of(analyte_2)) |>
        transmute(sample_id, {{ ratio_name }} := prot_1_symbol - prot_2_symbol)

      # Combine with metadata
      protein_meta_ratio_df <- inner_join(protein_ratio_metadata, ratio_df, by = "sample_id")

      # Filter out NA or Inf ratio values and check for what is left
      protein_meta_ratio_df <- protein_meta_ratio_df |>
        filter(!is.na(.data[[ratio_name]]), is.finite(.data[[ratio_name]]))
      if (nrow(protein_meta_ratio_df) < 1) {
        # Too few data points to fit a logistic model safely
        return(empty_result)
      }

      # Build formula and fit logistic regression model
      predictors <- c(ratio_name, other_predictors)
      formula <- reformulate(predictors, response = "pd")

      model <- glm(formula, data = protein_meta_ratio_df, family = binomial(link = "logit"))
      model_summary <- summary(model)

      tibble::tibble(
        protein_pair = protein_pair,
        mapped_analyte_1 = analyte_1,
        mapped_analyte_2 = analyte_2,
        prot_1_uniprot = prot_1_uniprot,
        prot_2_uniprot = prot_2_uniprot,
        gene_symbol1 = prot_1_symbol,
        gene_symbol2 = prot_2_symbol,
        estimate = model_summary$coefficients[ratio_name, "Estimate"],
        or_lower_ci = suppressMessages(stats::confint(model)[ratio_name, 1]),
        or_upper_ci = suppressMessages(stats::confint(model)[ratio_name, 2]),
        p_value = model_summary$coefficients[ratio_name, "Pr(>|z|)"],
        n_obs = nrow(model$model),
        converged = model$converged
      )
    },
    error = function(e) {
      empty_result
    }
  )

  return(result)
}
