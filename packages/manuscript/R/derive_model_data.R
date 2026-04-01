## SPDX-FileCopyrightText: 2026 Alkahest, Inc. <chronospd_data_access@alkahest.com>
##
## SPDX-License-Identifier: MIT

#' derive_model_data
#'
#' Derive model data by extracting relevant covariates and checking their presence
#' in metadata. Supports random-effects specified as **lmer-style character strings**.
#'
#' @param protein_id Character. The response column in `three_tibble_list$prot_data`.
#' @param three_tibble_list A list with components:
#'   \itemize{
#'     \item prot_data: tibble/data.frame with protein columns (incl. `protein_id`).
#'     \item meta_data: tibble/data.frame with covariates (same row order as prot_data).
#'     \item prot_annot: tibble/data.frame with annotation (used elsewhere).
#'   }
#' @param fixed_effects A list or character vector of fixed-effect RHS terms
#'   (e.g., `c("age", "sex", "x*y")`). Default empty means intercept-only.
#' @param random_effects Either:
#'   \itemize{
#'     \item **character** lmer-style RHS (e.g., `"1 | Donor + (0 + slope | Donor) + 1 | Plate"`), or
#'     \item **legacy list-of-lists** (only for `nlme` path).
#'   }
#'
#' @return A tibble containing the response column and all referenced covariates.
#' @export
derive_model_data <- function(protein_id, three_tibble_list, fixed_effects = list(), random_effects = NULL) {
  check_three_tibble_list(three_tibble_list)

  prot_data <- three_tibble_list$prot_data
  meta_data <- three_tibble_list$meta_data

  if (!protein_id %in% colnames(prot_data)) {
    stop(sprintf("protein_id '%s' not found in three_tibble_list$prot_data", protein_id))
  }

  cov_fixed <- extract_covariates_from_fixed(fixed_effects)
  cov_random <- if (is.character(random_effects)) extract_covariates_from_random_lmer(random_effects) else {
    if (is.list(random_effects) && length(random_effects) > 0) {
      vars <- c(
        unlist(lapply(random_effects, function(re) re$variables), use.names = FALSE),
        unlist(lapply(random_effects, function(re) strsplit(re$group, "[*:/|]")[[1]]), use.names = FALSE)
      )
      setdiff(unique(vars), c("0","1"))
    } else character()
  }

  covariates <- setdiff(unique(c(cov_fixed, cov_random)), c("0","1"))

  missing_cov <- setdiff(covariates, colnames(meta_data))
  if (length(missing_cov) > 0) {
    stop(sprintf("Some covariates are not present in meta_data: %s",
                 paste(missing_cov, collapse = ", ")))
  }

  model_data <- dplyr::bind_cols(
    prot_data |> dplyr::select(dplyr::all_of(protein_id)),
    meta_data |> dplyr::select(dplyr::all_of(covariates))
  )

  if (!is.data.frame(model_data)) {
    stop("model_data is not a data frame.")
  }
  model_data
}
