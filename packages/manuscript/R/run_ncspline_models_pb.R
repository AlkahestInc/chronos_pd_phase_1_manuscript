## SPDX-FileCopyrightText: 2026 Alkahest, Inc. <chronospd_data_access@alkahest.com>
##
## SPDX-License-Identifier: MIT

# natural cubic spline mixed model ####
# years_since_pd_onset + group + ns(years_since_pd_onset, df = DF):group
# Both groups: share one straight line
# PD:          straight line + nonlinear deviation
# Control:     just the straight line
run_ncspline_models_pb <- function(data, proteins, lme_results, df_range, mc.cores) {
  # Requires: pbmcapply, lmerTest, splines, dplyr

  # Function to process each protein
  process_protein_nc <- function(protein, data, lme_results, df_range) {
    # CRITICAL: Load packages in correct order in each worker to avoid ABI mismatch
    if (!isNamespaceLoaded("Matrix")) library(Matrix, quietly = TRUE)
    if (!isNamespaceLoaded("lme4")) library(lme4, quietly = TRUE)
    if (!isNamespaceLoaded("lmerTest")) library(lmerTest, quietly = TRUE)
    if (!isNamespaceLoaded("splines")) library(splines, quietly = TRUE)
    if (!require("dplyr", quietly = TRUE)) library(dplyr, quietly = TRUE)

    # Explicitly reference splines::ns to avoid namespace issues
    ns <- splines::ns

    # ---- locate protein entry in lme_results ----
    protein_idx <- which(sapply(lme_results, function(x) x$best_model_summary$protein$ID == protein))
    if (length(protein_idx) == 0) {
      return(NULL)
    }

    protein_info <- lme_results[[protein_idx]]$best_model_summary
    selected_formula <- protein_info$formula # original best model formula (string)
    lme_BIC_value <- protein_info$BIC # baseline LME BIC

    # ---- build base formula string with this protein as response (as in your code) ----
    base_formula <- as.formula(selected_formula)
    formula_str_base <- paste(deparse(base_formula), collapse = " ")
    # Set LHS to backticked protein index column (your original convention)
    formula_str_base <- sub("^[^~]+", paste0("`", protein, "`"), formula_str_base)

    # ---- table to collect per-df metrics ----
    protein_results <- data.frame(
      df = numeric(),
      n_knot = numeric(),
      LogLikelihood = numeric(),
      AIC = numeric(),
      BIC = numeric(),
      lme_BIC = numeric(),
      diff_lmeBIC = numeric(),
      max_vel_point = character(),
      stringsAsFactors = FALSE
    )

    # Keep track of best model only (avoids serialization issues with pbmclapply)
    best_model <- NULL
    best_BIC <- Inf
    best_df_val <- NA

    # ---- loop over candidate dfs (ML fits for fair comparison) ----
    for (i in seq_along(df_range)) {
      df_candidate <- df_range[i]

      # Replace ONLY the specific term:
      #   years_since_pd_onset * group
      # with:
      #   years_since_pd_onset + group + ns(years_since_pd_onset, df=DF):group
      new_term <- paste0("years_since_pd_onset + group + splines::ns(years_since_pd_onset, df=", df_candidate, "):group")
      new_formula_str <- gsub(
        pattern = "years_since_pd_onset\\s*\\*\\s*group",
        replacement = new_term,
        x = formula_str_base
      )
      new_formula <- as.formula(new_formula_str)

      model_ml <- tryCatch(
        lme4::lmer(new_formula, data = data, REML = FALSE),
        error = function(e) NULL
      )

      if (!is.null(model_ml)) {
        model_BIC <- BIC(model_ml)
        model_AIC <- AIC(model_ml)
        model_LL <- as.numeric(logLik(model_ml))
        current_diff_lmeBIC <- lme_BIC_value - model_BIC

        protein_results <- rbind(
          protein_results,
          data.frame(
            df = df_candidate,
            n_knot = df_candidate - 1,
            LogLikelihood = model_LL,
            AIC = model_AIC,
            BIC = model_BIC,
            lme_BIC = lme_BIC_value,
            diff_lmeBIC = current_diff_lmeBIC,
            max_vel_point = NA_character_,
            stringsAsFactors = FALSE
          )
        )

        # Keep track of best model for refitting
        if (model_BIC < best_BIC) {
          best_BIC <- model_BIC
          best_df_val <- df_candidate
          best_model <- model_ml
        }
      }
    }

    # If nothing fit, skip this protein
    if (nrow(protein_results) == 0 || is.null(best_model)) {
      return(NULL)
    }

    # ---- pick best df by smallest BIC ----
    best_row <- which.min(protein_results$BIC)
    best_ml <- best_model

    # ---- refit best model with REML for prediction ----
    best_reml <- tryCatch(
      lme4::lmer(update.formula(formula(best_ml), . ~ .), data = data, REML = TRUE),
      error = function(e) NULL
    )

    # ---- compute PD velocity (only for the best df row) ----
    if (!is.null(best_reml)) {
      # Prediction grid (same as your original; adjust as needed)
      new_data_nc <- expand.grid(
        years_since_pd_onset = seq(
          min(data$years_since_pd_onset, na.rm = TRUE),
          max(data$years_since_pd_onset, na.rm = TRUE),
          length.out = 200
        ),
        group = c("control", "pd")
      )
      new_data_nc <- new_data_nc %>%
        dplyr::mutate(
          gender = "male",
          age_at_baseline = median(data$age_at_baseline, na.rm = TRUE),
          recent_plasma_donations = "at least one donation",
          cluster = "1",
          race = "Caucasian"
        )
      new_data_nc$group <- factor(new_data_nc$group, levels = c("control", "pd"))

      # Marginal predictions (fixed effects only)
      new_data_nc$predicted <- predict(best_reml, newdata = new_data_nc, re.form = NA)

      new_data_pd <- subset(new_data_nc, group == "pd")
      new_data_control <- subset(new_data_nc, group == "control")

      max_vel_point <- NA_character_

      if (nrow(new_data_pd) >= 2) {
        pd_vel <- diff(new_data_pd$predicted) / diff(new_data_pd$years_since_pd_onset)
        # safe index (avoid off-by-one issues)
        pd_idx <- which.max(abs(pd_vel))
        pd_idx <- if (length(pd_idx) && pd_idx < nrow(new_data_pd)) pd_idx + 1L else nrow(new_data_pd)
        max_pd_value <- new_data_pd$years_since_pd_onset[pd_idx]
      } else {
        pd_vel <- numeric(0)
        max_pd_value <- NA_real_
      }

      if (nrow(new_data_control) >= 2) {
        control_vel <- diff(new_data_control$predicted) / diff(new_data_control$years_since_pd_onset)
        control_mag <- if (length(control_vel)) max(abs(control_vel), na.rm = TRUE) else -Inf
      } else {
        control_mag <- -Inf
      }

      pd_mag <- if (length(pd_vel)) max(abs(pd_vel), na.rm = TRUE) else -Inf

      if (!is.na(max_pd_value) && (is.infinite(control_mag) || pd_mag > control_mag)) {
        max_vel_point <- as.character(round(max_pd_value, 2))
      } else {
        max_vel_point <- "NA control"
      }

      # write back only into the best row
      protein_results$max_vel_point[best_row] <- max_vel_point
    }

    # Mark best row (optional, handy downstream)
    protein_results$is_best <- FALSE
    protein_results$is_best[best_row] <- TRUE

    protein_results <- protein_results |>
      dplyr::filter(is_best == TRUE) |>
      dplyr::select(-is_best)

    # return as before (ID is the numeric index key you used)
    return(list(ID = protein, results = protein_results))
  }

  # Run parallel processing with proper error handling
  results_list <- pbmcapply::pbmclapply(
    proteins,
    process_protein_nc,
    data = data,
    lme_results = lme_results,
    df_range = df_range,
    mc.cores = mc.cores,
    mc.preschedule = FALSE,
    mc.cleanup = TRUE,
    ignore.interactive = TRUE
  )

  # drop NULLs and name output by protein index as before
  results_list <- results_list[!sapply(results_list, is.null)]
  final_results_output <- setNames(
    lapply(results_list, function(x) x$results),
    sapply(results_list, function(x) x$ID)
  )
  return(final_results_output)
}
