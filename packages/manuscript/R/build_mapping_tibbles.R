## SPDX-FileCopyrightText: 2026 Alkahest, Inc. <chronospd_data_access@alkahest.com>
##
## SPDX-License-Identifier: MIT

#' Build analyte mapping tibbles from protein-pair tables
#'
#' Given a list of protein-pair tibbles, extracts analyte-to-UniProt
#' and gene symbol mappings for both proteins in each pair.
#'
#' @param protein_pair_tibbles List of tibbles containing protein pair data.
#' @param protein1_col Character. Column name for first protein ID.
#' @param protein2_col Character. Column name for second protein ID.
#' @param uniprot1_col Character. UniProt column for first protein.
#' @param uniprot2_col Character. UniProt column for second protein.
#' @param gene_symbol1_col Character. Gene symbol column for first protein.
#' @param gene_symbol2_col Character. Gene symbol column for second protein.
#'
#' @return A list of tibbles with standardized columns:
#'   platform_id, uniprot_id, gene_symbol
build_mapping_tibbles <- function(
  protein_pair_tibbles,
  protein1_col = "protein1",
  protein2_col = "protein2",
  uniprot1_col = "uniprot_id1",
  uniprot2_col = "uniprot_id2",
  gene_symbol1_col = "gene_symbol1",
  gene_symbol2_col = "gene_symbol2"
) {
  purrr::map(
    protein_pair_tibbles,
    function(tbl) {
      tbl1 <- tbl |>
        dplyr::select(
          platform_id = dplyr::all_of(protein1_col),
          uniprot_id = dplyr::all_of(uniprot1_col),
          gene_symbol = dplyr::all_of(gene_symbol1_col)
        )

      tbl2 <- tbl |>
        dplyr::select(
          platform_id = dplyr::all_of(protein2_col),
          uniprot_id = dplyr::all_of(uniprot2_col),
          gene_symbol = dplyr::all_of(gene_symbol2_col)
        )

      list(tbl1, tbl2)
    }
  ) |>
    purrr::flatten()
}
