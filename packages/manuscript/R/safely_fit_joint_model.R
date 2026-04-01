## SPDX-FileCopyrightText: 2026 Alkahest, Inc. <chronospd_data_access@alkahest.com>
##
## SPDX-License-Identifier: MIT

#' Safely fit a single joint model
#' Force random slope! JMbayes2 requires random slopes when single longitudinal outcome
#' @keywords internal
safely_fit_joint_model <- function(data, spec, protein) {
  tryCatch(
    {
      cli::cli_alert_info(glue::glue("Preparing data for protein {protein}"))

      long_vars <- all.vars(spec$longitudinal_formula)
      sel_vars <- c(data$id_var, data$time_var, protein, long_vars)

      protein_data <- data$longitudinal |>
        dplyr::select(dplyr::all_of(sel_vars)) |>
        na.omit()

      cli::cli_alert_info(glue::glue(
        "Protein {protein}: n_long = {nrow(protein_data)}, ",
        "n_id = {dplyr::n_distinct(protein_data[[data$id_var]])}"
      ))

      cli::cli_alert_info(glue::glue("Protein {protein}: fitting lme at {Sys.time()}"))
      lme_fit <- nlme::lme(
        fixed   = update(spec$longitudinal_formula, paste(protein, "~ .")),
        random  = as.formula(spec$random_formula),
        data    = protein_data,
        control = nlme::lmeControl(opt = "optim", maxIter = 200, returnObject = TRUE)
      )
      cli::cli_alert_success(glue::glue("Protein {protein}: finished lme at {Sys.time()}"))

      cli::cli_alert_info(glue::glue("Protein {protein}: fitting Cox model at {Sys.time()}"))
      cox_fit <- survival::coxph(
        formula = spec$survival_formula,
        data    = data$survival
      )
      cli::cli_alert_success(glue::glue("Protein {protein}: finished Cox model at {Sys.time()}"))

      functional_form <- build_functional_form(protein, spec$association)

      cli::cli_alert_info(glue::glue("Protein {protein}: starting JMbayes2::jm at {Sys.time()}"))
      jm_fit <- JMbayes2::jm(
        Surv_object      = cox_fit,
        Mixed_objects    = lme_fit,
        time_var         = data$time_var,
        functional_forms = functional_form,
        n_iter           = spec$mcmc_settings$iter,
        n_burnin         = spec$mcmc_settings$warmup,
        n_thin           = spec$mcmc_settings$thin,
        n_chains         = spec$mcmc_settings$chains
      )
      cli::cli_alert_success(glue::glue("Protein {protein}: finished JMbayes2::jm at {Sys.time()}"))

      list(
        protein     = protein,
        model       = jm_fit,
        status      = "success",
        convergence = check_convergence(jm_fit),
        summary     = extract_model_summary(jm_fit, protein)
      )
    },
    error = function(e) {
      cli::cli_alert_danger(
        glue::glue("JM failed for protein {protein}: {conditionMessage(e)}")
      )
      list(
        protein     = protein,
        model       = NULL,
        status      = "failed",
        error       = as.character(e),
        convergence = NULL,
        summary     = NULL
      )
    }
  )
}
