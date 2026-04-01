## SPDX-FileCopyrightText: 2026 Alkahest, Inc. <chronospd_data_access@alkahest.com>
##
## SPDX-License-Identifier: MIT

#' Create long format data for forest plots
#'
#' @param protein_pairs Vector of protein pair IDs to process
#' @param chronos_full Chronos data with columns: protein_pair, cohort_label, window_start, window_end, odds_ratio, or_lower_ci, or_upper_ci, p_value
#' @param ukb_filtered UKB data with columns: protein_pair, cohort_label, window_start, window_end, odds_ratio, or_lower_ci, or_upper_ci, p_value
#' @param ppmi_full PPMI data with columns: protein_pair, estimate, odds_ratio, or_lower_ci, or_upper_ci, p_value
#' @param gnpc_cohortC_full GNPC Cohort C data with columns: protein_pair, estimate, odds_ratio, or_lower_ci, or_upper_ci, p_value
#' @param gnpc_cohortJ_full GNPC Cohort J data with columns: protein_pair, estimate, odds_ratio, or_lower_ci, or_upper_ci, p_value
#' @param gnpc_cohortL_full GNPC Cohort L data with columns: protein_pair, estimate, odds_ratio, or_lower_ci, or_upper_ci, p_value
#'
#' @return Data frame with forest plot data for all protein pairs across all cohorts
create_forest_data <- function(protein_pairs, chronos_full, ukb_filtered,
                                ppmi_full, gnpc_cohortC_full,
                                gnpc_cohortJ_full, gnpc_cohortL_full) {
  forest_data_list <- list()

  for (pair in protein_pairs) {
    # Chronos data for this pair
    chronos_forest <- chronos_full |>
      dplyr::filter(protein_pair == pair) |>
      dplyr::arrange(desc(window_start)) |>
      dplyr::mutate(
        cohort = glue::glue("Chronos_{window_start}_{window_end}"),
        cohort_group = "Chronos"
      ) |>
      dplyr::select(
        protein_pair, cohort, cohort_group, window_start, window_end,
        odds_ratio, or_lower_ci, or_upper_ci, p_value
      )

    # UKB data for this pair (filtered windows, ordered most recent first for forest plot)
    ukb_forest <- ukb_filtered |>
      dplyr::filter(protein_pair == pair) |>
      dplyr::mutate(
        cohort = glue::glue("UKB_{window_start}_{window_end}"),
        cohort_group = "UKB"
      ) |>
      dplyr::arrange(desc(window_start)) |>
      dplyr::select(
        protein_pair, cohort, cohort_group, window_start, window_end,
        odds_ratio, or_lower_ci, or_upper_ci, p_value
      )

    # Ensure both UKB windows are present in correct order
    expected_ukb <- data.frame(
      protein_pair = pair,
      cohort = c("UKB_-5_-1", "UKB_-9_-5"),
      cohort_group = "UKB",
      window_start = c(-5, -9),
      window_end = c(-1, -5)
    )
    ukb_forest <- expected_ukb |>
      dplyr::left_join(ukb_forest |> dplyr::select(-protein_pair, -cohort_group, -window_start, -window_end),
                       by = "cohort"
      )

    # PPMI data for this pair
    ppmi_forest <- ppmi_full |>
      dplyr::filter(protein_pair == pair) |>
      dplyr::mutate(
        cohort = "PPMI",
        cohort_group = "PPMI",
        window_start = NA,
        window_end = NA
      ) |>
      dplyr::select(
        protein_pair, cohort, cohort_group, window_start, window_end,
        odds_ratio, or_lower_ci, or_upper_ci, p_value
      )

    # If no PPMI data, create empty row
    if (nrow(ppmi_forest) == 0) {
      ppmi_forest <- data.frame(
        protein_pair = pair,
        cohort = "PPMI",
        cohort_group = "PPMI",
        window_start = NA,
        window_end = NA,
        odds_ratio = NA,
        or_lower_ci = NA,
        or_upper_ci = NA,
        p_value = NA
      )
    }

    # GNPC cohorts for this pair (keep most significant)
    gnpc_c_forest <- gnpc_cohortC_full |>
      dplyr::filter(protein_pair == pair) |>
      keep_most_significant() |>
      dplyr::mutate(
        cohort = "GNPC C",
        cohort_group = "GNPC",
        window_start = NA,
        window_end = NA
      ) |>
      dplyr::select(
        protein_pair, cohort, cohort_group, window_start, window_end,
        odds_ratio, or_lower_ci, or_upper_ci, p_value
      )

    # If no GNPC C data, create empty row
    if (nrow(gnpc_c_forest) == 0) {
      gnpc_c_forest <- data.frame(
        protein_pair = pair,
        cohort = "GNPC C",
        cohort_group = "GNPC",
        window_start = NA,
        window_end = NA,
        odds_ratio = NA,
        or_lower_ci = NA,
        or_upper_ci = NA,
        p_value = NA
      )
    }

    gnpc_j_forest <- gnpc_cohortJ_full |>
      dplyr::filter(protein_pair == pair) |>
      keep_most_significant() |>
      dplyr::mutate(
        cohort = "GNPC J",
        cohort_group = "GNPC",
        window_start = NA,
        window_end = NA
      ) |>
      dplyr::select(
        protein_pair, cohort, cohort_group, window_start, window_end,
        odds_ratio, or_lower_ci, or_upper_ci, p_value
      )

    # If no GNPC J data, create empty row
    if (nrow(gnpc_j_forest) == 0) {
      gnpc_j_forest <- data.frame(
        protein_pair = pair,
        cohort = "GNPC J",
        cohort_group = "GNPC",
        window_start = NA,
        window_end = NA,
        odds_ratio = NA,
        or_lower_ci = NA,
        or_upper_ci = NA,
        p_value = NA
      )
    }

    gnpc_l_forest <- gnpc_cohortL_full |>
      dplyr::filter(protein_pair == pair) |>
      keep_most_significant() |>
      dplyr::mutate(
        cohort = "GNPC L",
        cohort_group = "GNPC",
        window_start = NA,
        window_end = NA
      ) |>
      dplyr::select(
        protein_pair, cohort, cohort_group, window_start, window_end,
        odds_ratio, or_lower_ci, or_upper_ci, p_value
      )

    # If no GNPC L data, create empty row
    if (nrow(gnpc_l_forest) == 0) {
      gnpc_l_forest <- data.frame(
        protein_pair = pair,
        cohort = "GNPC L",
        cohort_group = "GNPC",
        window_start = NA,
        window_end = NA,
        odds_ratio = NA,
        or_lower_ci = NA,
        or_upper_ci = NA,
        p_value = NA
      )
    }

    # Combine all cohorts for this pair (reverse order for coord_flip display)
    pair_data <- dplyr::bind_rows(
      gnpc_l_forest,
      gnpc_j_forest,
      gnpc_c_forest,
      ppmi_forest,
      ukb_forest,
      chronos_forest
    )

    forest_data_list[[pair]] <- pair_data
  }

  # Combine all pairs
  dplyr::bind_rows(forest_data_list)
}
