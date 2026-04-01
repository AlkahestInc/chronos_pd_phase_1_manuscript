## SPDX-FileCopyrightText: 2026 Alkahest, Inc. <chronospd_data_access@alkahest.com>
##
## SPDX-License-Identifier: MIT

#' Analyze multiple protein ratios from a protein-pair tibble
#'
#' Iterates over protein pairs, expands analyte combinations,
#' fits logistic regression models for each ratio, and returns
#' a combined results table with adjusted p-values.
#'
#' @param protein_pair_tibble Tibble with protein pair definitions.
#' @param protein_ratio_metadata Tibble with sample metadata.
#' @param protein_intensity_log2 Tibble of log2 protein intensities.
#' @param other_predictors Character vector of covariate column names.
#'
#' @return Tibble of model results with ORs and BH-adjusted p-values.
analyze_multiple_ratio_by_tibble <- function(
  protein_pair_tibble,
  protein_ratio_metadata,
  protein_intensity_log2,
  other_predictors = c("age_at_visit", "sex")
) {
  results <- purrr::pmap_dfr(
    protein_pair_tibble,
    function(protein1, protein2, protein_pair, ...) {
      analytes_1 <- process_protein_ratio_prob_str(
        probe = protein1,
        mapping_cross_olink_platforms = TRUE,
        mapping_cross_somalogic_platforms = TRUE
      )
      analytes_2 <- process_protein_ratio_prob_str(protein2,
        mapping_cross_olink_platforms = TRUE,
        mapping_cross_somalogic_platforms = TRUE
      )
      prot_1_uniprot <- paste(unique(get_uniprot_by_prob(protein1,
        mapping_cross_olink_platforms = TRUE,
        mapping_cross_somalogic_platforms = TRUE
      )), sep = ":")
      prot_2_uniprot <- paste(unique(get_uniprot_by_prob(protein2,
        mapping_cross_olink_platforms = TRUE,
        mapping_cross_somalogic_platforms = TRUE
      )), sep = ":")
      prot_1_symbol <- paste(unique(get_gene_symbol_by_prob(protein1,
        mapping_cross_olink_platforms = TRUE,
        mapping_cross_somalogic_platforms = TRUE
      )), sep = ":")
      prot_2_symbol <- paste(unique(get_gene_symbol_by_prob(protein2,
        mapping_cross_olink_platforms = TRUE,
        mapping_cross_somalogic_platforms = TRUE
      )), sep = ":")

      analyte_pairs <- expand.grid(
        analyte_1 = analytes_1,
        analyte_2 = analytes_2,
        stringsAsFactors = FALSE
      )

      purrr::pmap_dfr(
        analyte_pairs,
        ~ analyze_single_ratio(
          ..1,
          ..2,
          protein_pair,
          prot_1_uniprot,
          prot_2_uniprot,
          prot_1_symbol,
          prot_2_symbol,
          protein_ratio_metadata,
          protein_intensity_log2,
          other_predictors
        )
      )
    }
  )

  # Add additional calculated columns
  results <- results |>
    dplyr::mutate(or = exp(estimate)) |>
    dplyr::arrange(p_value)

  return(results)
}
