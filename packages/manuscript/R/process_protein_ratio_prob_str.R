## SPDX-FileCopyrightText: 2026 Alkahest, Inc. <chronospd_data_access@alkahest.com>
##
## SPDX-License-Identifier: MIT

#' Process protein ratio probe string
#'
#' Normalizes Somalogic probe IDs or, for non-Somalogic probes,
#' retrieves annotation column names via UniProt mapping.
#'
#' @param probe Character. Analyte or probe ID.
#' @param mapping_cross_olink_platforms Bool. Indicate if to map across
#'        different sets of Olink analytes. If FALSE, use Olink 5k analytes.
#'        Default to TRUE
#' @param mapping_cross_somalogic_platforms Bool. Indicate if to map across
#'        different sets of somalogic analytes. If FALSE, use Soma 7k analytes.
#'        Default to TRUE
#'
#' @return Character vector of processed probe strings or annotation column names.
process_protein_ratio_prob_str <- function(
  probe,
  mapping_cross_olink_platforms = TRUE,
  mapping_cross_somalogic_platforms = TRUE
) {
  if (grepl("^seq", probe)) {
    gsub("\\.", "_", probe)
  }
  uniprot_id <- get_uniprot_by_prob(probe,
    mapping_cross_olink_platforms = mapping_cross_olink_platforms,
    mapping_cross_somalogic_platforms = mapping_cross_somalogic_platforms
  )
  annotations_by_prot <- annotations |>
    dplyr::filter(uni_prot %in% uniprot_id) |>
    dplyr::pull(column_name)
}
