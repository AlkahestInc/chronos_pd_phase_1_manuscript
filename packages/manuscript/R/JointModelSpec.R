## SPDX-FileCopyrightText: 2026 Alkahest, Inc. <chronospd_data_access@alkahest.com>
##
## SPDX-License-Identifier: MIT

#' R6 Class for Joint Model Specification
#'
#' Modern OOP approach for model configuration
#' @export
JointModelSpec <- R6::R6Class(
  "JointModelSpec",
  public = list(
    # Public fields
    longitudinal_formula = NULL,
    survival_formula = NULL,
    random_formula = NULL,
    association = "value",

    # Constructor
    initialize = function(longitudinal_formula = NULL,
                          survival_formula = NULL,
                          random_formula = ~ 1 | id,
                          association = c("value", "slope", "area", "shared_random")) {
      self$longitudinal_formula <- longitudinal_formula
      self$survival_formula <- survival_formula
      self$random_formula <- random_formula
      self$association <- match.arg(association, several.ok = TRUE)

      private$validate()
    },

    # Modern method chaining
    set_longitudinal = function(formula) {
      self$longitudinal_formula <- formula
      private$validate()
      invisible(self)
    },
    set_survival = function(formula) {
      self$survival_formula <- formula
      private$validate()
      invisible(self)
    },
    set_association = function(type) {
      self$association <- match.arg(type, c("value", "slope", "area", "shared_random"),
                                    several.ok = TRUE
      )
      invisible(self)
    },

    # MCMC configuration with modern defaults
    set_mcmc = function(chains = 4, iter = 5000, warmup = 1000,
                        thin = 2, adapt_delta = 0.95) {
      private$mcmc_config <- list(
        chains = chains,
        iter = iter,
        warmup = warmup,
        thin = thin,
        adapt_delta = adapt_delta
      )
      invisible(self)
    },
    print = function() {
      cli::cli_h1("Joint Model Specification")

      cli::cli_h2("Formulas")
      if (!is.null(self$longitudinal_formula)) {
        cli::cli_li("Longitudinal: {.code {deparse(self$longitudinal_formula)}}")
      }
      if (!is.null(self$survival_formula)) {
        cli::cli_li("Survival: {.code {deparse(self$survival_formula)}}")
      }
      cli::cli_li("Random: {.code {deparse(self$random_formula)}}")

      cli::cli_h2("Association")
      cli::cli_li("Type(s): {.val {self$association}}")

      cli::cli_h2("MCMC Settings")
      mcmc <- private$mcmc_config
      cli::cli_li("{.val {mcmc$chains}} chains × {.val {mcmc$iter}} iterations")
      cli::cli_li("Warmup: {.val {mcmc$warmup}}, Thin: {.val {mcmc$thin}}")

      invisible(self)
    }
  ),
  private = list(
    mcmc_config = list(
      chains = 4,
      iter = 5000,
      warmup = 1000,
      thin = 2,
      adapt_delta = 0.95
    ),
    validate = function() {
      if (!is.null(self$longitudinal_formula) && !inherits(self$longitudinal_formula, "formula")) {
        cli::cli_abort("Longitudinal formula must be a formula object")
      }
      if (!is.null(self$survival_formula) && !inherits(self$survival_formula, "formula")) {
        cli::cli_abort("Survival formula must be a formula object")
      }
    }
  ),
  active = list(
    # Computed properties
    is_complete = function() {
      !is.null(self$longitudinal_formula) && !is.null(self$survival_formula)
    },
    mcmc_settings = function() {
      private$mcmc_config
    }
  )
)
