## SPDX-FileCopyrightText: 2026 Alkahest, Inc. <chronospd_data_access@alkahest.com>
##
## SPDX-License-Identifier: MIT

here::i_am("src/40_manuscript/main_figures/figure_2/figure_2j.R")

gene_symbols_to_plot <- c("DDC", "SUMF1", "ITGAV", "CD276", "QPCT", "OIP5", "TUBA1A")

chronos_dataset_files <- list(
  "SomaScan" = glue::glue("{get_config()$input_dir}/name=chronos_somalogic_all-pd_master_summary/version=29806/dataset.parquet"),
  "Explore HT" = glue::glue("{get_config()$input_dir}/name=chronos_olink_all-pd_master_summary/version=29809/dataset.parquet"),
  "NULISA CNS" = glue::glue("{get_config()$input_dir}/name=chronos_alamar_all-pd_master_summary/version=29812/dataset.parquet")
)

gnpc_dataset_files <- list(
  "GNPC C" = glue::glue("{get_config()$input_dir}/name=GNPC_C_single_protein_validation_result/version=none/dataset.parquet"),
  "GNPC J" = glue::glue("{get_config()$input_dir}/name=GNPC_J_single_protein_validation_result/version=none/dataset.parquet"),
  "GNPC L" = glue::glue("{get_config()$input_dir}/name=GNPC_L_single_protein_validation_result/version=none/dataset.parquet")
)

ukb_dataset_files <- list(
  "UKB" = glue::glue("{get_config()$input_dir}/name=ukb_single_protein_adjusted_linear_regression/version=20251216_060400_+0000/dataset.parquet")
)

ppmi_dataset_files <- list(
  "PPMI" = glue::glue("{get_config()$input_dir}/name=PPMI_pd_hits_lme_results/version=Feb8_2026/dataset.parquet")
)

olink_assay_id2_gene_symbol_file <- glue::glue("{get_config()$input_dir}/name=protein_platform_uniprot_ids/version=29404/dataset.parquet")

output_figure_file <- glue::glue("{get_config()$output_dir}/main_figures/figure_2j.pdf")
dir.create(dirname(output_figure_file), showWarnings = FALSE, recursive = TRUE)


olink_assay_id2_gene_symbol <- olink_assay_id2_gene_symbol_file |>
  arrow::open_dataset() |>
  dplyr::collect() |>
  tidyr::separate_longer_delim(gene_symbol, delim = ":") |>
  dplyr::filter(platform == "olink") |>
  dplyr::distinct(
    gene_symbol,
    uniprot_id,
  )

#####
# Load CHRONOS data
#####
chronos_data <- purrr::map(
  chronos_dataset_files,
  ~ {
    .x |>
      arrow::open_dataset() |>
      dplyr::select(
        probe_id,
        uniprot_id,
        gene_symbol,
        estimate,
        std_err,
        cohens_d,
        p_val,
        q_val
      ) |>
      dplyr::collect()
  }
) |>
  dplyr::bind_rows(.id = "platform") |>
  dplyr::mutate(
    phase_sample = "all",
    cohort = "Chronos-PD"
  ) |>
  tidyr::separate_longer_delim(uniprot_id, delim = "|") |>
  tidyr::separate_longer_delim(gene_symbol, delim = "|") |>
  dplyr::group_by(uniprot_id) |>
  dplyr::slice_min(p_val, n = 1) |>
  dplyr::ungroup()

#####
# Load GNPC data
#####
gnpc_data <- purrr::map(
  gnpc_dataset_files,
  ~ {
    .x |>
      arrow::open_dataset() |>
      dplyr::select(
        uniprot_id = prot_uniprot,
        gene_symbol,
        estimate = pdYes_estimate,
        p_val = pdYes_p_value,
        std_err = pdYes_std_error
      ) |>
      dplyr::collect() |>
      dplyr::group_by(uniprot_id) |>
      dplyr::slice_min(p_val, n = 1) |>
      dplyr::ungroup()
  }
) |>
  dplyr::bind_rows(.id = "cohort") |>
  dplyr::distinct()

#####
# Load UKB data
#####
ukb_data <- purrr::map(
  ukb_dataset_files,
  ~ {
    .x |>
      arrow::open_dataset() |>
      dplyr::select(
        gene_symbol = Protein,
        estimate = beta,
        p_val = p.value,
        std_err = std.error
      ) |>
      dplyr::collect() |>
      tidyr::separate_longer_delim(gene_symbol, delim = "_")
  }
) |>
  dplyr::bind_rows(.id = "cohort") |>
  #####
  # Add gene_symbol and restrict data to those we have uniprot_id <-> gene_symbol mappings
  #####
  dplyr::inner_join(
    olink_assay_id2_gene_symbol,
    by = dplyr::join_by(gene_symbol)
  ) |>
  dplyr::select(
    cohort, uniprot_id, gene_symbol, estimate, p_val, std_err
  ) |>
  dplyr::group_by(uniprot_id) |>
  dplyr::slice_min(p_val, n = 1) |>
  dplyr::ungroup() |>
  dplyr::distinct()

#####
# Load PPMI data
#####
ppmi_data <- purrr::map(
  ppmi_dataset_files,
  ~ {
    .x |>
      arrow::open_dataset() |>
      dplyr::select(
        uniprot_id = UNIPROT,
        gene_symbol = SYMBOL,
        estimate = Value,
        p_val = p.value,
        std_err = Std.Error
      ) |>
      dplyr::collect()
  }
) |>
  dplyr::bind_rows(.id = "cohort") |>
  tidyr::separate_longer_delim(gene_symbol, delim = "|") |>
  dplyr::group_by(uniprot_id) |>
  dplyr::slice_min(p_val, n = 1) |>
  dplyr::ungroup() |>
  dplyr::distinct()

master_summary <- dplyr::bind_rows(
  chronos_data |> dplyr::select(cohort, uniprot_id, gene_symbol, estimate, p_val, std_err),
  ppmi_data,
  gnpc_data,
  ukb_data
) |>
  dplyr::filter(
    gene_symbol %in% gene_symbols_to_plot
  ) |>
  dplyr::rename(
    log2FC = estimate
  ) |>
  dplyr::mutate(
    lower_ci = log2FC - 1.96 * std_err,
    upper_ci = log2FC + 1.96 * std_err,
    gene_symbol = factor(gene_symbol, levels = gene_symbols_to_plot)
  )

p <- master_summary |>
  dplyr::mutate(
    colour = pval_to_stars(p_val)
  ) |>
  # Only plot the most significant gene_symbol per cohort
  dplyr::group_by(cohort, gene_symbol) |>
  dplyr::slice_min(p_val) |>
  dplyr::ungroup() |>
  protein_forest_plot(
    cohort_order = c("GNPC L", "GNPC J", "GNPC C", "PPMI", "UKB", "Chronos-PD")
  )

ggplot2::ggsave(
  filename = output_figure_file,
  plot = p,
  width = 12.5, height = 3.0, units = "in",
  dpi = 600,
  bg = "white"
)
