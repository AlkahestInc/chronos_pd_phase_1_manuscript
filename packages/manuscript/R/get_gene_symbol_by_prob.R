## SPDX-FileCopyrightText: 2026 Alkahest, Inc. <chronospd_data_access@alkahest.com>
##
## SPDX-License-Identifier: MIT

#' Get gene symbol from analyte/probe ID
#'
#' Given a probe/analyte ID, look up the corresponding gene symbol
#' using the appropriate platform-specific mapping table.
#'
#' @param probe Character. Analyte or probe ID.
#' @param mapping_cross_olink_platforms Bool. Indicate if to map across
#'        different sets of Olink analytes. If FALSE, use Olink 5k analytes.
#'        Default to TRUE
#' @param mapping_cross_somalogic_platforms Bool. Indicate if to map across
#'        different sets of somalogic analytes. If FALSE, use Soma 7k analytes.
#'        Default to TRUE
#'
#' @return Character vector of gene symbols.
get_gene_symbol_by_prob <- function(
  probe,
  mapping_cross_olink_platforms = TRUE,
  mapping_cross_somalogic_platforms = TRUE,
  platform_mapping_file = "/home/rstudio/input_data/metadata/protein_platform_uniprot_ids.csv"
) {
  # Load platform mapping file
  platform_mapping <- readr::read_csv(platform_mapping_file)

  olink_analyte_mapping <- if (mapping_cross_olink_platforms) {
    platform_mapping |>
      dplyr::filter(platform == "olink")
  } else {
    platform_mapping |>
      dplyr::filter(platform == "olink", platform_version == "5k")
  }

  somalogic_analyte_mapping <- if (mapping_cross_somalogic_platforms) {
    platform_mapping |>
      dplyr::filter(platform == "somalogic")
  } else {
    platform_mapping |>
      dplyr::filter(platform == "somalogic", platform_version == "7k")
  }

  if (grepl("^OID", probe)) {
    olink_analyte_mapping |>
      dplyr::filter(
        platform_id == probe
      ) |>
      dplyr::pull(gene_symbol)
  } else if (grepl("^seq", probe)) {
    somalogic_analyte_mapping |>
      dplyr::filter(
        platform_id == probe
      ) |>
      dplyr::pull(gene_symbol)

  } else {
    other_analyte_mapping |>
      dplyr::filter(
        platform_id == probe
      ) |>
      dplyr::pull(gene_symbol)
  }
}
