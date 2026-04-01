## SPDX-FileCopyrightText: 2026 Alkahest, Inc. <chronospd_data_access@alkahest.com>
##
## SPDX-License-Identifier: MIT

#' Identify significant diagnosis codes for a given stage
#'
#' @param res_chronos Chronos results data frame
#' @param res_us US results data frame
#' @param stage Stage to filter ("before" or "after")
#' @param max_codes Maximum number of codes to return (default: NULL for no limit / all codes)
#'
#' @return A data frame with significant codes and their plot order
identify_significant_codes_for_stage <- function(res_chronos, res_us, stage, max_codes = NULL) {
  # Identify codes to plot - significant codes in Chronos-PD
  significant_chronos <- res_chronos |>
    dplyr::filter(FDR < 0.05) |>
    dplyr::filter(stage == !!stage) |>
    dplyr::arrange(dplyr::desc(logOR)) |>
    dplyr::mutate(order_to_plot = dplyr::row_number()) |>
    dplyr::select(DIAGNOSIS_CODE_CATEGORY, order_to_plot)

  # Retrieve the ORs for those same codes in US Established-PD
  significant_us <- res_us |>
    dplyr::filter(FDR < 0.05) |>
    dplyr::filter(stage == !!stage) |>
    dplyr::filter(DIAGNOSIS_CODE_CATEGORY %in%
                    significant_chronos$DIAGNOSIS_CODE_CATEGORY) |>
    dplyr::arrange(dplyr::desc(logOR))

  if (!is.null(max_codes)) {
    significant_us <- significant_us |>
      utils::head(max_codes)
  }

  significant_us <- significant_us |>
    dplyr::mutate(order_to_plot = dplyr::row_number()) |>
    dplyr::select(DIAGNOSIS_CODE_CATEGORY, order_to_plot)

  return(significant_us)
}
