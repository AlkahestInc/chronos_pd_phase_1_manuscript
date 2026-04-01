## SPDX-FileCopyrightText: 2026 Alkahest, Inc. <chronospd_data_access@alkahest.com>
##
## SPDX-License-Identifier: MIT

#' Keep most significant result per protein pair
#'
#' @param data Data frame containing protein pair results
#' @param protein_pair_col Column name for protein pairs (default: "protein_pair")
#' @param p_value_col Column name for p-values (default: "p_value")
#'
#' @return Data frame with only the most significant result per protein pair
keep_most_significant <- function(data, protein_pair_col = "protein_pair", p_value_col = "p_value") {
  data |>
    dplyr::group_by(dplyr::across(dplyr::all_of(protein_pair_col))) |>
    dplyr::slice_min(order_by = .data[[p_value_col]], n = 1, with_ties = FALSE) |>
    dplyr::ungroup()
}
