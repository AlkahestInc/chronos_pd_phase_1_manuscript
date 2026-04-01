## SPDX-FileCopyrightText: 2026 Alkahest, Inc. <chronospd_data_access@alkahest.com>
##
## SPDX-License-Identifier: MIT

#' Load and combine PD data from Chronos and US cohorts
#'
#' @param chronos_res Chronos-PD Phase 1 results
#' @param us_res US Established-PD results
#'
#' @return A data frame with combined results from both cohorts
load_and_combine_odds_ratio_data <- function(chronos_res, us_res) {

  # Select relevant columns
  chronos_res <- chronos_res |>
    dplyr::select(
      stage,
      DIAGNOSIS_CODE_CATEGORY,
      A,
      B,
      C,
      D,
      OR,
      SE = StdError,
      PValue,
      FDR,
      logOR,
      LowerCI,
      UpperCI,
      significant,
      label
    )

  us_res <- us_res |>
    dplyr::select(
      stage,
      DIAGNOSIS_CODE_CATEGORY,
      A,
      B,
      C,
      D,
      OR,
      SE = StdError,
      PValue,
      FDR,
      logOR,
      LowerCI,
      UpperCI,
      significant,
      label
    )

  # Combine into one data.frame
  all_res <- dplyr::bind_rows(chronos_res, us_res,
                               .id = "dataset"
  ) |>
    dplyr::mutate(dataset = dplyr::case_when(
      dataset == 1 ~ "Chronos-PD",
      dataset == 2 ~ "US Population-Level PD"
    ))

  # Filter out Z codes for visualization purposes
  all_res <- all_res[!grepl("^Z", all_res$DIAGNOSIS_CODE_CATEGORY), ]

  return(list(
    all_res = all_res,
    chronos_res = chronos_res,
    us_res = us_res
  ))
}
