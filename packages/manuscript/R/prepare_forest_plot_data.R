## SPDX-FileCopyrightText: 2026 Alkahest, Inc. <chronospd_data_access@alkahest.com>
##
## SPDX-License-Identifier: MIT

#' Prepare forest plot data for a given stage
#'
#' @param all_res Combined results data frame
#' @param significant_codes Data frame with significant codes and plot order
#' @param stage Stage to filter ("before" or "after")
#'
#' @return A data frame ready for plotting
prepare_forest_plot_data <- function(all_res, significant_codes, stage) {
  stage_label <- dplyr::case_when(
    stage == "after" ~ "Post-diagnosis",
    stage == "before" ~ "Pre-diagnosis"
  )

  plotting_data <- all_res |>
    dplyr::mutate(label = stringr::str_to_sentence(label)) |>
    dplyr::mutate(stage2 = dplyr::case_when(
      stage == "after" ~ "Post-diagnosis",
      stage == "before" ~ "Pre-diagnosis"
    )) |>
    dplyr::filter(stage2 == stage_label) |>
    dplyr::mutate(stage2 = factor(stage2, levels = c("Pre-diagnosis", "Post-diagnosis"))) |>
    dplyr::inner_join(significant_codes, by = "DIAGNOSIS_CODE_CATEGORY") |>
    dplyr::mutate(label_to_plot = paste0(DIAGNOSIS_CODE_CATEGORY, ": ", label))

  return(plotting_data)
}
