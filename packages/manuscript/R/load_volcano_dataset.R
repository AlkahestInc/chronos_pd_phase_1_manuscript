## SPDX-FileCopyrightText: 2026 Alkahest, Inc. <chronospd_data_access@alkahest.com>
##
## SPDX-License-Identifier: MIT

#' Load and process differential expression dataset for volcano plot
#'
#' @param input_dataset_file Path to the input parquet dataset
#' @param paper_pd_list A data frame containing PD hit proteins with uniprot_id column
#'
#' @return A processed data frame with gene_symbol, uniprot_id, log2_fc_all, p_val_all, q_val_all, and point_group
load_volcano_dataset <- function(input_dataset_file, paper_pd_list) {
  dataset_result <- arrow::open_dataset(
    input_dataset_file
  ) |>
    dplyr::select(
      gene_symbol,
      uniprot_id,
      log2_fc_all,
      p_val_all,
      q_val_all
    ) |>
    dplyr::distinct() |>
    dplyr::mutate(
      point_group = dplyr::case_when(
        uniprot_id %in% paper_pd_list[["uniprot_id"]] & p_val_all < 0.05 ~ "PD_hit",
        q_val_all < 0.05 ~ "FDR_sig",
        p_val_all < 0.05 ~ "P_sig",
        TRUE ~ "NS"
      )
    ) |>
    dplyr::collect()

  return(dataset_result)
}
