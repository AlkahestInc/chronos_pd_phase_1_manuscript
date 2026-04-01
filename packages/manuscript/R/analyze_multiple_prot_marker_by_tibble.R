## SPDX-FileCopyrightText: 2026 Alkahest, Inc. <chronospd_data_access@alkahest.com>
##
## SPDX-License-Identifier: MIT

#' Analyze multiple single protein markers from a tibble with above single
#' protein analysis function
#'
#' Iterates over protein probes, resolves analytes, fits regression models,
#' optionally adjusts p-values, and returns a sorted result table.
#'
#' @param single_prot_tibble Tibble. Contains platform_id column.
#' @param single_prot_metadata Tibble. Sample metadata.
#' @param protein_intensity_log2 Tibble. Log2 protein intensities.
#' @param annotations Tibble. Annotation table.
#' @param calc_adj_p Logical. Whether to compute BH-adjusted p-values.
#' @param sort Logical. Whether to sort by R-squared.
#' @param other_predictors Character vector. Covariates for regression.
#' @param mapping_cross_olink_platforms Bool. Indicate if to map across
#'        different sets of Olink analytes during the mapping step.
#'        If FALSE, use Soma 7k analytes. Default to FALSE
#' @param mapping_cross_somalogic_platforms Bool. Indicate if to map across
#'        different sets of somalogic analytes during the mapping step.
#'        If FALSE, use Soma 7k analytes. Default to FALSE
#'
#' @return A tibble of regression results for all probes.
analyze_multiple_prot_marker_by_tibble <- function(
  single_prot_tibble,
  single_prot_metadata,
  protein_intensity_log2,
  annotations,
  calc_adj_p = TRUE,
  sort = TRUE,
  mapping_cross_olink_platforms = FALSE,
  mapping_cross_somalogic_platforms = FALSE
) {
  other_predictors <- c("age_at_visit", "sex", "pd")

  results <- purrr::pmap_dfr(
    list(single_prot_tibble$platform_id),
    function(probe_id) {
      analytes <- process_protein_marker_prob_str(probe_id,
        mapping_cross_olink_platforms = mapping_cross_olink_platforms,
        mapping_cross_somalogic_platforms = mapping_cross_somalogic_platforms
      )
      if (is.vector(analytes) && length(analytes) > 1) {
        analytes <- as.list(analytes)
      } else if (!(is.character(analytes) && length(analytes) == 1)) {
        analytes <- as.list("")
      } else {
        analytes <- as.list(analytes)
      }

      prot_uniprot <- paste(unique(get_uniprot_by_prob(probe_id,
        mapping_cross_olink_platforms = mapping_cross_olink_platforms,
        mapping_cross_somalogic_platforms = mapping_cross_somalogic_platforms
      )), sep = ":")
      gene_symbol <- paste(unique(get_gene_symbol_by_prob(probe_id,
        mapping_cross_olink_platforms = mapping_cross_olink_platforms,
        mapping_cross_somalogic_platforms = mapping_cross_somalogic_platforms
      )), sep = ":")

      purrr::map_dfr(
        analytes,
        ~ analyze_single_protein_marker(
          probe_id,
          ..1,
          prot_uniprot,
          gene_symbol,
          single_prot_metadata,
          protein_intensity_log2,
          annotations,
          other_predictors
        )
      )
    }
  )

  if (calc_adj_p) {
    # BH-adjusted p-values per term
    results <- results |>
      dplyr::mutate(
        dplyr::across(
          tidyselect::ends_with("_p_value"),
          ~ stats::p.adjust(.x, method = "BH"),
          .names = "{.col}_bh_adj"
        )
    )
  }
  if (sort) {
    results <- results |>
      dplyr::arrange(dplyr::desc(r_squared))
  }

  # ---- reorder columns ----
  model_cols <- c("protein", "r_squared", "adj_r_squared", "n")
  cn <- colnames(results)
  pred_cols <- setdiff(cn, model_cols)
  # all other predictor columns sorted alphabetically
  other_cols <- sort(pred_cols[!grepl("^protein_", pred_cols)])
  final_cols <- c(model_cols, other_cols)

  results <- results |> dplyr::select(tidyselect::all_of(final_cols))

  return(results)
}
