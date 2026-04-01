## SPDX-FileCopyrightText: 2026 Alkahest, Inc. <chronospd_data_access@alkahest.com>
##
## SPDX-License-Identifier: MIT

run_linear_spline_models_pb <- function(data, proteins, lme_results, knot_range, mc.cores) {

  # Function to process each protein
  process_protein <- function(protein, data, lme_results, knot_range) {
    # CRITICAL: Load packages in correct order in each worker to avoid ABI mismatch
    if (!isNamespaceLoaded("Matrix")) library(Matrix, quietly = TRUE)
    if (!isNamespaceLoaded("lme4")) library(lme4, quietly = TRUE)
    if (!isNamespaceLoaded("lmerTest")) library(lmerTest, quietly = TRUE)
    if (!isNamespaceLoaded("lspline")) library(lspline, quietly = TRUE)

    idx <- which(sapply(lme_results, function(x) x$best_model_summary$protein$ID == protein))
    if (length(idx) == 0) {
      return(list(ID = protein, results = NULL))
    } # Skip if protein not found

    selected_formula <- lme_results[[idx]]$best_model_summary$formula
    lme_BIC_value <- lme_results[[idx]]$best_model_summary$BIC

    # Convert the selected formula to a formula object and deparse to a string.
    base_formula <- as.formula(selected_formula)
    formula_str_base <- paste(deparse(base_formula), collapse = " ")
    # Wrap protein in backticks so that names with dashes are valid variable names.
    formula_str_base <- sub("^[^~]+", paste0("`", protein, "`"), formula_str_base)

    # Data frame to store results for each candidate knot for this protein.
    protein_results <- data.frame(
      Knot = numeric(),
      AIC_ls = numeric(),
      BIC_ls = numeric(),
      BIC_lme = numeric(),
      diff_lsBIC = numeric(),
      group_diff_sig = character(),
      stringsAsFactors = FALSE
    )

    # Loop over candidate knots.
    for (knot in knot_range) {
      ## 1. Linear Spline Model (using lspline) with REML = FALSE for AIC/BIC.
      new_term_lin <- paste0("lspline(years_since_pd_onset, ", knot, ") * group")
      new_formula_str_lin <- gsub(
        "years_since_pd_onset\\s*\\*\\s*group",
        new_term_lin, formula_str_base
      )
      formula_lin <- as.formula(new_formula_str_lin)

      model_lin <- tryCatch(
        lmerTest::lmer(formula_lin, data = data, REML = FALSE),
        error = function(e) {
          message("Error fitting model for ", protein, " knot ", knot, ": ", e$message)
          return(NULL)
        }
      )

      if (!is.null(model_lin)) {
        AIC_ls <- AIC(model_lin)
        BIC_ls <- BIC(model_lin)
        # Evaluate group difference significance using a separate REML=TRUE fit.
        model_lin_REML <- tryCatch(
          lmerTest::lmer(formula_lin, data = data, REML = TRUE),
          error = function(e) {
            return(NULL)
          }
        )

        if (!is.null(model_lin_REML)) {
          coeff_table <- as.data.frame(summary(model_lin_REML)$coefficients)
          pattern1 <- paste0("lspline\\(years_since_pd_onset,\\s*", knot, "\\)1:group")
          pattern2 <- paste0("lspline\\(years_since_pd_onset,\\s*", knot, "\\)2:group")
          pval1 <- min(coeff_table[grep(pattern1, rownames(coeff_table)), "Pr(>|t|)"], na.rm = TRUE)
          pval2 <- min(coeff_table[grep(pattern2, rownames(coeff_table)), "Pr(>|t|)"], na.rm = TRUE)
          if (!is.na(pval1) & !is.na(pval2)) {
            if (pval1 < 0.05 & pval2 < 0.05) {
              group_diff_sig <- "both"
            } else if (pval1 < 0.05 & pval2 >= 0.05) {
              group_diff_sig <- "first"
            } else if (pval1 >= 0.05 & pval2 < 0.05) {
              group_diff_sig <- "second"
            } else {
              group_diff_sig <- "neither"
            }
          } else {
            group_diff_sig <- NA
          }
        } else {
          group_diff_sig <- NA
        }
      } else {
        AIC_ls <- NA
        BIC_ls <- NA
        group_diff_sig <- NA
      }

      ## 2. Compute differences in BIC.
      diff_lsBIC <- lme_BIC_value - BIC_ls

      # Append candidate knot results for this protein.
      protein_results <- rbind(protein_results, data.frame(
        Knot = knot,
        AIC_ls = AIC_ls,
        BIC_ls = BIC_ls,
        BIC_lme = lme_BIC_value,
        diff_lsBIC = diff_lsBIC,
        group_diff_sig = group_diff_sig,
        stringsAsFactors = FALSE
      ))
    }

    return(list(ID = protein, results = protein_results))
  }

  # Run parallel processing with proper error handling
  results <- pbmcapply::pbmclapply(
    proteins,
    process_protein,
    data = data,
    lme_results = lme_results,
    knot_range = knot_range,
    mc.cores = mc.cores,
    mc.preschedule = FALSE,
    mc.cleanup = TRUE
  )

  # Filter out NULL results
  results <- Filter(Negate(is.null), results)

  final_results <- setNames(
    lapply(results, function(x) x$results),
    sapply(results, function(x) x$ID)
  )

  return(final_results)
}
