## SPDX-FileCopyrightText: 2026 Alkahest, Inc. <chronospd_data_access@alkahest.com>
##
## SPDX-License-Identifier: MIT

#' compare_multiple_models
#'
#' Compare multiple combinations of fixed/random effects and choose the best by BIC.
#' **Random effects should be provided as lmer-style RHS strings** (character vector).
#'
#' @param protein_id Character.
#' @param three_tibble_list See `derive_model_data()`.
#' @param fixed_effects_list A character vector of formula-like RHS strings
#'   (e.g., `c("x*y + sex + age", "x*y + sex")`) **or** a list where each element is a
#'   character vector of tokens to be collapsed with `" + "`. Each element corresponds to one model.
#' @param random_effects_list Character vector of lmer-style random RHS strings,
#'   e.g., `c("1 | Donor + 1 | Plate", "1 + slope | Donor + 1 | Plate", ...)`.
#' @param package "lmerTest" (default) or "lmer". (Use legacy lists for nlme.)
#' @param mc.cores Integer. If >1 and non-Windows, uses parallel::mclapply.
#' @param ... Passed to `explain_by_LME`.
#'
#' @return A list with: models (BICs), best_model, best_model_fixed, best_model_random, best_bic
#' @export
compare_multiple_models <- function(
    protein_id,
    three_tibble_list,
    fixed_effects_list,
    random_effects_list,
    package = "lmerTest",
    mc.cores = 1,
    ...
) {
  if (!is.character(random_effects_list)) {
    stop("random_effects_list must be a character vector of lmer-style RHS strings.")
  }
  if (!package %in% c("lmer", "lmerTest")) {
    stop("For random-effects strings, 'package' must be 'lmer' or 'lmerTest'. For 'nlme' use legacy list-of-lists.")
  }

  fixed_vec <- normalize_fixed_list(fixed_effects_list)

  grid <- tidyr::expand_grid(
    fixed = fixed_vec,
    random = random_effects_list
  ) |>
    dplyr::mutate(model_id = dplyr::row_number())

  run_one <- function(row) {
    fe <- row$fixed
    re <- row$random
    fit <- tryCatch(
      explain_by_LME(protein_id, three_tibble_list,
                     fixed_effects = fe,
                     random_effects = re,
                     package = package, ...),
      error = function(e) NULL
    )
    if (is.null(fit)) return(list(BIC = NA_real_))
    fit_ml <- tryCatch(refit_mixed(fit, REML = FALSE), error = function(e) fit)
    list(BIC = stats::BIC(fit_ml))
  }

  apply_fun <- if (mc.cores > 1 && .Platform$OS.type != "windows") {
    function(idx, FUN) parallel::mclapply(idx, FUN, mc.cores = mc.cores)
  } else {
    function(idx, FUN) lapply(idx, FUN)
  }

  res_list <- apply_fun(seq_len(nrow(grid)), function(i) run_one(grid[i, , drop = FALSE]))
  BICs <- vapply(res_list, function(x) x$BIC, numeric(1))

  grid$BIC <- BICs
  grid <- grid |> dplyr::arrange(BIC)

  best <- grid |> dplyr::slice(1)

  best_model <- tryCatch(
    explain_by_LME(
      protein_id, three_tibble_list,
      fixed_effects = best$fixed,    # RHS string directly
      random_effects = best$random,
      package = package,
      ...
    ),
    error = function(e) NULL
  )

  list(
    models = purrr::imap(
      split(grid, seq_len(nrow(grid))),
      ~ list(BIC = .x$BIC, fixed = .x$fixed, random = .x$random)
    ),
    best_model = best_model,
    best_model_fixed = best$fixed,      # keep RHS string
    best_model_random = best$random,
    best_bic = best$BIC
  )
}
