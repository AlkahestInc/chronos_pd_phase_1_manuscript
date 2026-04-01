## SPDX-FileCopyrightText: 2026 Alkahest, Inc. <chronospd_data_access@alkahest.com>
##
## SPDX-License-Identifier: MIT

#' Run logistic regression for a single protein
#'
#' @param gene_symbol Character string. Gene symbol for the protein.
#' @param analysis_data Data frame containing demographic and clinical variables.
#' @param analysis_prot Data frame containing protein expression data with columns as GNPC analytes.
#'
#' @return A tibble with one row containing:
#'   - protein_name: Protein identifier
#'   - gene_symbol: Gene symbol
#'   - estimate: Beta coefficient for the protein
#'   - or_lower_ci: Lower bound of 95% confidence interval for odds ratio
#'   - or_upper_ci: Upper bound of 95% confidence interval for odds ratio
#'   - p_value: P-value for the protein coefficient
#'   - n_obs: Number of observations used in the model
#'   - converged: Logical indicating if model converged
#'
#' @examples
#' \dontrun{
#' result <- analyze_single_protein("NTRK3", analysis_data, analysis_prot)
#' }
analyze_single_protein <- function(analyte, protein_name, gene_symbol, analysis_data, analysis_prot) {
  empty_results <- tibble::tibble(
    analyte = analyte,
    protein_name = protein_name,
    gene_symbol = gene_symbol,
    estimate = NA_real_,
    or_lower_ci = NA_real_,
    or_upper_ci = NA_real_,
    p_value = NA_real_,
    n_obs = NA_integer_,
    converged = FALSE
  )

  # Construct column name
  col_name <- protein_name

  # Check if column exists
  if (!col_name %in% names(analysis_prot)) {
    return(empty_results)
  }

  # Create a copy of the data to avoid modifying original, and add protein values
  selected_prot_tibble <- analysis_prot[c("sample_id", col_name)] |> rename(!!gene_symbol := col_name)
  model_data <- inner_join(analysis_data, selected_prot_tibble, by = "sample_id")

  # Check if protein has valid values
  if (all(is.na(model_data[[gene_symbol]]))) {
    return(empty_results)
  }

  # Ensure pd is properly factored
  model_data$pd <- factor(model_data$pd,
    levels = c("No", "Yes")
  )

  # Build formula dynamically
  formula_str <- paste0("pd ~ ", gene_symbol, " + age_at_visit + sex")
  model_formula <- stats::as.formula(formula_str)

  # Try to fit the model with error handling
  tryCatch(
    {
      model <- stats::glm(model_formula, data = model_data, family = "binomial")

      # Check convergence
      converged <- model$converged

      # Extract model summary
      model_summary <- summary(model)

      # Get coefficient information for the protein
      estimate_val <- model_summary$coefficients[gene_symbol, "Estimate"]
      pval <- model_summary$coefficients[gene_symbol, "Pr(>|z|)"]

      # Calculate confidence intervals
      conf_int <- suppressMessages(stats::confint(model)[gene_symbol, ])

      # Return results
      tibble::tibble(
        analyte = analyte,
        protein_name = protein_name,
        gene_symbol = gene_symbol,
        estimate = estimate_val,
        or_lower_ci = conf_int[1],
        or_upper_ci = conf_int[2],
        p_value = pval,
        n_obs = nrow(model$model),
        converged = converged
      )
    },
    error = function(e) {
      empty_results
    }
  )
}
