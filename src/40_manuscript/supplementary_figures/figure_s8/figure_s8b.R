## SPDX-FileCopyrightText: 2026 Alkahest, Inc. <chronospd_data_access@alkahest.com>
##
## SPDX-License-Identifier: MIT

here::i_am("src/40_manuscript/supplementary_figures/figure_s8/figure_s8b.R")

gene_symbol <- "CXCL12"

input_proteomics_files <- list(
  somalogic = glue::glue("{get_config()$input_dir}/name=Chronos_phase1_proteomics_somalogic_analysis_ready_three_tibble/version=29650/dataset.rds"),
  olink = glue::glue("{get_config()$input_dir}/name=Chronos_phase1_proteomics_olink_analysis_ready_three_tibble/version=29465/dataset.rds"),
  alamar = glue::glue("{get_config()$input_dir}/name=Chronos_phase1_proteomics_alamar_analysis_ready_three_tibble/version=29649/dataset.rds"),
  biognosys = glue::glue("{get_config()$input_dir}/name=Chronos_phase1_proteomics_biognosys_analysis_ready_three_tibble/version=29639/dataset.rds")
)

platform_2_uniprot_mapping_file <- glue::glue("{get_config()$input_dir}/name=protein_platform_uniprot_ids/version=29404/dataset.parquet")

output_figure_file <- glue::glue("{get_config()$output_dir}/supplementary_figures/figure_s8/figure_s8b.pdf")
dir.create(dirname(output_figure_file), showWarnings = FALSE, recursive = TRUE)

platform_2_uniprot_mapping <- platform_2_uniprot_mapping_file |>
  arrow::open_dataset() |>
  dplyr::filter(
    gene_symbol %in% !!gene_symbol
  ) |>
  dplyr::select(
    gene_symbol,
    probe_id = platform_id,
    uniprot_id
  ) |>
  dplyr::arrange(probe_id) |>
  dplyr::distinct() |>
  dplyr::collect()

uniprot_id <- platform_2_uniprot_mapping |>
  dplyr::distinct(uniprot_id) |>
  dplyr::pull()

.biognosys_gene_symbol2label <- function(gene_symbol) {
  biognosys_lookup_file <- glue::glue("{get_config()$input_dir}/name=chronos_biognosys_all-pd_master_summary/version=29864/dataset.parquet")

  biognosys_lookup <- biognosys_lookup_file |>
    arrow::open_dataset() |>
    dplyr::select(gene_symbol, probe_id, label) |>
    dplyr::filter(
      gene_symbol %in% !!gene_symbol
    ) |>
    dplyr::collect()
  biognosys_rename_vector <- setNames(biognosys_lookup$label, biognosys_lookup$probe_id)

  return(biognosys_rename_vector)
}

target_columns <- c(.biognosys_gene_symbol2label(gene_symbol), setNames(platform_2_uniprot_mapping[["probe_id"]], platform_2_uniprot_mapping[["probe_id"]]))

############################################
#### Load/Process Three-Tibble Datasets ####
############################################
message("Loading proteomics datasets")
proteomics_data <- purrr::map(
  input_proteomics_files,
  ~ {
    three_tbl <- readr::read_rds(.x)

    # Find matching columns
    matching_probes <- intersect(
      colnames(three_tbl$prot_data),
      target_columns
    )

    # Return NULL if the file has no matches for the requested probes
    if (length(matching_probes) == 0) {
      return(NULL)
    }

    tbl <- three_tbl$meta_data |>
      # Get the minimum metadata required for this script
      dplyr::select(
        ID,
        sample_id
      ) |>
      dplyr::inner_join(
        three_tbl$prot_data |>
          # Only get a subset of data for the proteins we need in this dataset
          dplyr::select(
            ID,
            dplyr::any_of(target_columns)
          ),
        by = "ID"
      ) |>
      dplyr::select(-ID) |>
      tidyr::pivot_longer(
        cols = -sample_id,
        names_to = "probe_id",
      )

    return(tbl)
  },
  progress = TRUE
) |>
  dplyr::bind_rows(.id = "platform") |>
  dplyr::left_join(
    platform_2_uniprot_mapping,
    by = "probe_id"
  )

plotting_data <- proteomics_data |>
  tidyr::pivot_wider(
    id_cols = sample_id,
    names_from = probe_id,
    values_from = value
  ) |>
  dplyr::select(-sample_id)

p <- plotting_data |>
  GGally::ggpairs(
    lower = list(continuous = lower_smooth),
    upper = list(continuous = upper_cor_spearman),
    diag = list(continuous = GGally::wrap("barDiag", bins = 30)),
    title = glue::glue("{gene_symbol} - {uniprot_id}")
  ) +
  ggplot2::theme_minimal()


ggplot2::ggsave(
  filename = output_figure_file,
  plot = p,
  width = 8.0,
  height = 8.0,
  units = "in",
  dpi = 600,
  bg = "white"
)
