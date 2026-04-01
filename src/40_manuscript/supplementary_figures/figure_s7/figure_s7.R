## SPDX-FileCopyrightText: 2026 Alkahest, Inc. <chronospd_data_access@alkahest.com>
##
## SPDX-License-Identifier: MIT

here::i_am("src/40_manuscript/supplementary_figures/figure_s7/figure_s7.R")

# Define input files
input_files <- list(
  chronos_list1 = glue::glue("{get_config()$input_dir}/name=chronos_list1/version=202510/dataset.parquet"),
  ppmi_list1 = glue::glue("{get_config()$input_dir}/name=ppmi_list1/version=202510/dataset.parquet"),
  ukb_ratio = glue::glue("{get_config()$input_dir}/name=ukb_ratio/version=20251208_205951_+0000/dataset.parquet"),
  ukb_single = glue::glue("{get_config()$input_dir}/name=ukb_single/version=20260207_075402_+0000/dataset.parquet"),
  gnpc_c_ratio = glue::glue("{get_config()$input_dir}/name=gnpc_c_ratio/version=202510/dataset.parquet"),
  gnpc_c_single = glue::glue("{get_config()$input_dir}/name=gnpc_c_single/version=202510/dataset.parquet"),
  gnpc_j_ratio = glue::glue("{get_config()$input_dir}/name=gnpc_j_ratio/version=202510/dataset.parquet"),
  gnpc_j_single = glue::glue("{get_config()$input_dir}/name=gnpc_j_single/version=202510/dataset.parquet"),
  gnpc_l_ratio = glue::glue("{get_config()$input_dir}/name=gnpc_l_ratio/version=202510/dataset.parquet"),
  gnpc_l_single = glue::glue("{get_config()$input_dir}/name=gnpc_l_single/version=202510/dataset.parquet")
)
assay_id2uniprot_id_file <- glue::glue("{get_config()$input_dir}/name=protein_platform_uniprot_ids/version=29404/dataset.parquet")

# Define output file
output_figure_file <- glue::glue("{get_config()$output_dir}/supplementary_figures/figure_s7/figure_s7.pdf")
dir.create(dirname(output_figure_file), showWarnings = FALSE, recursive = TRUE)

pair_id_order <- c(
  "OID44508_o_seq.9191.8",
  "seq.31524.5_o_seq.9191.8",
  "OID43875_o_seq.8275.31",
  "OID45216_o_OID45275",
  "OID45275_o_seq.20579.50",
  "OID45195_o_OID45275",
  "seq.20093.9_o_seq.9191.8",
  "seq.3805.16_o_seq.9191.8",
  "seq.10746.24_o_seq.9191.8",
  "seq.5694.57_o_seq.9191.8"
)

cohort_order <- c(
  "GNPC C",
  "GNPC J",
  "GNPC L",
  "PPMI",
  "UKB [-5,-1)",
  "UKB [-9,-5)",
  "Chronos [-1,1)",
  "Chronos [-2,0)",
  "Chronos [-3,-1)",
  "Chronos [-4,-2)",
  "Chronos [-5,-3)",
  "Chronos [-6,-4)",
  "Chronos [-7,-5)"
)

lookup <- assay_id2uniprot_id_file |>
  arrow::open_dataset() |>
  dplyr::distinct(
    platform_id,
    gene_symbol
  ) |>
  dplyr::collect()


pairs_to_plot <- tibble::tibble(
  pair_id = pair_id_order
) |>
  tidyr::separate(
    pair_id,
    into = c("id.numerator", "id.denominator"),
    sep = "_o_",
    remove = FALSE
  ) |>
  # Join gene_symbol for numerator
  dplyr::left_join(
    lookup,
    by = c("id.numerator" = "platform_id")
  ) |>
  dplyr::rename(gene_symbol.numerator = gene_symbol) |>
  # Join gene_symbol for denominator
  dplyr::left_join(
    lookup,
    by = c("id.denominator" = "platform_id")
  ) |>
  dplyr::rename(gene_symbol.denominator = gene_symbol) |>
  dplyr::mutate(
    facet_label = glue::glue("{gene_symbol.numerator} / {gene_symbol.denominator}")
  ) |>
  tidyr::pivot_longer(
    cols = c(
      id.numerator, id.denominator,
      gene_symbol.numerator, gene_symbol.denominator
    ),
    names_to = c(".value", "type"),
    names_sep = "\\."
  ) |>
  dplyr::select(pair_id, type, id, gene_symbol, facet_label) |>
  dplyr::mutate(
    pair_id = forcats::fct_relevel(pair_id, pair_id_order)
  ) |> dplyr::arrange(pair_id) |>
  dplyr::mutate(
    facet_label = forcats::fct_inorder(facet_label)
  )


#####
# Load the protein ratio statistical data
#####
ratio_data <- list()
ratio_data[["Chronos"]] <- input_files[["chronos_list1"]] |>
  arrow::open_dataset() |>
  dplyr::filter(
    analysis_type == "ratio"
  ) |>
  dplyr::collect() |>
  dplyr::mutate(
    x = odds_ratio,
    y = glue::glue("Chronos [{window_start},{window_end})"),
    ymin = or_lower_ci,
    ymax = or_upper_ci,
    p_value,
    pair_id = ratio,
    .keep = "none"
  ) |>
  dplyr::inner_join(
    pairs_to_plot,
    by = "pair_id"
  ) |>
  dplyr::mutate(
    type = "ratio"
  ) |>
  dplyr::distinct() |>
  dplyr::filter(!is.na(x)) |>
  dplyr::group_by(
    y, type, facet_label
  ) |>
  dplyr::slice_min(p_value, n = 1) |>
  dplyr::filter(
    y %in% cohort_order
  ) |>
  # Restrict to only those pairs we were initially interested in
  dplyr::semi_join(
    pairs_to_plot,
    by = "pair_id"
  )

ratio_data[["PPMI"]] <- input_files[["ppmi_list1"]] |>
  arrow::open_dataset() |>
  dplyr::filter(
    analysis_type == "ratio"
  ) |>
  dplyr::collect() |>
  dplyr::mutate(
    x = or,
    y = "PPMI",
    ymin = exp(or_lower_ci),
    ymax = exp(or_upper_ci),
    p_value,
    pair_id = protein_pair,
    .keep = "none"
  ) |>
  dplyr::inner_join(
    pairs_to_plot,
    by = "pair_id"
  ) |>
  dplyr::mutate(
    type = "ratio"
  ) |>
  dplyr::distinct() |>
  dplyr::filter(!is.na(x)) |>
  dplyr::group_by(
    y, type, facet_label
  ) |>
  dplyr::slice_min(p_value, n = 1) |>
  dplyr::filter(
    y %in% cohort_order
  ) |>
  # Restrict to only those pairs we were initially interested in
  dplyr::semi_join(
    pairs_to_plot,
    by = "pair_id"
  )

ratio_data[["GNPC C"]] <- input_files[["gnpc_c_ratio"]] |>
  arrow::open_dataset() |>
  dplyr::collect() |>
  dplyr::mutate(
    x = or,
    y = "GNPC C",
    ymin = exp(or_lower_ci),
    ymax = exp(or_upper_ci),
    p_value,
    pair_id = protein_pair,
    .keep = "none"
  ) |>
  dplyr::inner_join(
    pairs_to_plot,
    by = "pair_id"
  ) |>
  dplyr::mutate(
    type = "ratio"
  ) |>
  dplyr::distinct() |>
  dplyr::filter(!is.na(x)) |>
  dplyr::group_by(
    y, type, facet_label
  ) |>
  dplyr::slice_min(p_value, n = 1) |>
  dplyr::filter(
    y %in% cohort_order
  ) |>
  # Restrict to only those pairs we were initially interested in
  dplyr::semi_join(
    pairs_to_plot,
    by = "pair_id"
  )

ratio_data[["GNPC J"]] <- input_files[["gnpc_j_ratio"]] |>
  arrow::open_dataset() |>
  dplyr::collect() |>
  dplyr::mutate(
    x = or,
    y = "GNPC J",
    ymin = exp(or_lower_ci),
    ymax = exp(or_upper_ci),
    p_value,
    pair_id = protein_pair,
    .keep = "none"
  ) |>
  dplyr::inner_join(
    pairs_to_plot,
    by = "pair_id"
  ) |>
  dplyr::mutate(
    type = "ratio"
  ) |>
  dplyr::distinct() |>
  dplyr::filter(!is.na(x)) |>
  dplyr::group_by(
    y, type, facet_label
  ) |>
  dplyr::slice_min(p_value, n = 1) |>
  dplyr::filter(
    y %in% cohort_order
  ) |>
  # Restrict to only those pairs we were initially interested in
  dplyr::semi_join(
    pairs_to_plot,
    by = "pair_id"
  )

ratio_data[["GNPC L"]] <- input_files[["gnpc_l_ratio"]] |>
  arrow::open_dataset() |>
  dplyr::collect() |>
  dplyr::mutate(
    x = or,
    y = "GNPC L",
    ymin = exp(or_lower_ci),
    ymax = exp(or_upper_ci),
    p_value,
    pair_id = protein_pair,
    .keep = "none"
  ) |>
  dplyr::inner_join(
    pairs_to_plot,
    by = "pair_id"
  ) |>
  dplyr::mutate(
    type = "ratio"
  ) |>
  dplyr::distinct() |>
  dplyr::filter(!is.na(x)) |>
  dplyr::group_by(
    y, type, facet_label
  ) |>
  dplyr::slice_min(p_value, n = 1) |>
  dplyr::filter(
    y %in% cohort_order
  ) |>
  # Restrict to only those pairs we were initially interested in
  dplyr::semi_join(
    pairs_to_plot,
    by = "pair_id"
  )

ratio_data[["UKB"]] <- input_files[["ukb_ratio"]] |>
  arrow::open_dataset() |>
  dplyr::collect() |>
  dplyr::mutate(
    x = odds_ratio,
    y = glue::glue("UKB [{window_start},{window_end})"),
    ymin = conf.low,
    ymax = conf.high,
    p_value = p.value,
    pair_id,
    .keep = "none"
  ) |>
  dplyr::inner_join(
    pairs_to_plot,
    by = "pair_id"
  ) |>
  dplyr::mutate(
    type = "ratio"
  ) |>
  dplyr::distinct() |>
  dplyr::filter(!is.na(x)) |>
  dplyr::group_by(
    y, type, facet_label
  ) |>
  dplyr::slice_min(p_value, n = 1) |>
  dplyr::filter(
    y %in% cohort_order
  ) |>
  # Restrict to only those pairs we were initially interested in
  dplyr::semi_join(
    pairs_to_plot,
    by = "pair_id"
  )


#####
# Individual protein data
#####
individual_data <- list()
individual_data[["Chronos"]] <- input_files[["chronos_list1"]] |>
  arrow::open_dataset() |>
  dplyr::filter(
    analysis_type == "individual"
  ) |>
  dplyr::collect() |>
  dplyr::mutate(
    id = protein_name,
    x = odds_ratio,
    y = glue::glue("Chronos [{window_start},{window_end})"),
    ymin = or_lower_ci,
    ymax = or_upper_ci,
    p_value,
    .keep = "none"
  ) |>
  dplyr::inner_join(
    pairs_to_plot,
    by = "id"
  ) |>
  dplyr::distinct() |>
  dplyr::filter(!is.na(x)) |>
  dplyr::group_by(
    y, type, facet_label
  ) |>
  dplyr::slice_min(p_value, n = 1) |>
  dplyr::filter(
    y %in% cohort_order
  ) |>
  # Restrict to only those pairs we were initially interested in
  dplyr::semi_join(
    pairs_to_plot,
    by = "pair_id"
  )

individual_data[["PPMI"]] <- input_files[["ppmi_list1"]] |>
  arrow::open_dataset() |>
  dplyr::filter(
    analysis_type == "individual"
  ) |>
  dplyr::collect() |>
  dplyr::mutate(
    gene_symbol,
    x = or,
    y = "PPMI",
    ymin = exp(or_lower_ci),
    ymax = exp(or_upper_ci),
    p_value,
    .keep = "none"
  ) |>
  dplyr::inner_join(
    pairs_to_plot,
    by = "gene_symbol"
  ) |>
  dplyr::distinct() |>
  dplyr::filter(!is.na(x)) |>
  dplyr::group_by(
    y, type, facet_label
  ) |>
  dplyr::slice_min(p_value, n = 1) |>
  dplyr::filter(
    y %in% cohort_order
  ) |>
  # Restrict to only those pairs we were initially interested in
  dplyr::semi_join(
    pairs_to_plot,
    by = "pair_id"
  )

individual_data[["GNPC C"]] <- input_files[["gnpc_c_single"]] |>
  arrow::open_dataset() |>
  dplyr::collect() |>
  dplyr::mutate(
    gene_symbol,
    x = or,
    y = "GNPC C",
    ymin = exp(or_lower_ci),
    ymax = exp(or_upper_ci),
    p_value,
    .keep = "none"
  ) |>
  dplyr::inner_join(
    pairs_to_plot,
    by = "gene_symbol"
  ) |>
  dplyr::distinct() |>
  dplyr::filter(!is.na(x)) |>
  dplyr::group_by(
    y, type, facet_label
  ) |>
  dplyr::slice_min(p_value, n = 1) |>
  dplyr::filter(
    y %in% cohort_order
  ) |>
  # Restrict to only those pairs we were initially interested in
  dplyr::semi_join(
    pairs_to_plot,
    by = "pair_id"
  )

individual_data[["GNPC J"]] <- input_files[["gnpc_j_single"]] |>
  arrow::open_dataset() |>
  dplyr::collect() |>
  dplyr::mutate(
    gene_symbol,
    x = or,
    y = "GNPC J",
    ymin = exp(or_lower_ci),
    ymax = exp(or_upper_ci),
    p_value,
    .keep = "none"
  ) |>
  dplyr::inner_join(
    pairs_to_plot,
    by = "gene_symbol"
  ) |>
  dplyr::distinct() |>
  dplyr::filter(!is.na(x)) |>
  dplyr::group_by(
    y, type, facet_label
  ) |>
  dplyr::slice_min(p_value, n = 1) |>
  dplyr::filter(
    y %in% cohort_order
  ) |>
  # Restrict to only those pairs we were initially interested in
  dplyr::semi_join(
    pairs_to_plot,
    by = "pair_id"
  )

individual_data[["GNPC L"]] <- input_files[["gnpc_l_single"]] |>
  arrow::open_dataset() |>
  dplyr::collect() |>
  dplyr::mutate(
    gene_symbol,
    x = or,
    y = "GNPC L",
    ymin = exp(or_lower_ci),
    ymax = exp(or_upper_ci),
    p_value,
    .keep = "none"
  ) |>
  dplyr::inner_join(
    pairs_to_plot,
    by = "gene_symbol"
  ) |>
  dplyr::distinct() |>
  dplyr::filter(!is.na(x)) |>
  dplyr::group_by(
    y, type, facet_label
  ) |>
  dplyr::slice_min(p_value, n = 1) |>
  dplyr::filter(
    y %in% cohort_order
  ) |>
  # Restrict to only those pairs we were initially interested in
  dplyr::semi_join(
    pairs_to_plot,
    by = "pair_id"
  )

individual_data[["UKB"]] <- input_files[["ukb_single"]] |>
  arrow::open_dataset() |>
  dplyr::collect() |>
  dplyr::mutate(
    gene_symbol = Protein,
    x = odds_ratio,
    y = glue::glue("UKB [{window_start},{window_end})"),
    ymin = conf.low,
    ymax = conf.high,
    p_value = p.value,
    .keep = "none"
  ) |>
  dplyr::inner_join(
    pairs_to_plot,
    by = "gene_symbol"
  ) |>
  dplyr::distinct() |>
  dplyr::filter(!is.na(x)) |>
  dplyr::group_by(
    y, type, facet_label
  ) |>
  dplyr::slice_min(p_value, n = 1) |>
  dplyr::filter(
    y %in% cohort_order
  ) |>
  # Restrict to only those pairs we were initially interested in
  dplyr::semi_join(
    pairs_to_plot,
    by = "pair_id"
  )



.protein_forest_plot <- function(dataset) {
  y_levels <- levels(dataset$y)

  stripe_df <- dataset |>
    dplyr::distinct(cohort, y) |>
    dplyr::mutate(
      position = match(as.character(y), y_levels)
    ) |>
    dplyr::group_by(cohort) |>
    dplyr::summarise(
      xmin = min(position) - 0.5,
      xmax = max(position) + 0.5,
      .groups = "drop"
    ) |>
    dplyr::arrange(xmin) |>
    dplyr::mutate(
      fill = rep_len(c("#FFFFFF", "#DDDDDD"), dplyr::n())
    )

  p <- dataset |>
    ggplot2::ggplot(
      ggplot2::aes(
        x = y,
        y = x,
        colour = colour
      )
    ) +
    ggplot2::geom_rect(
      data = stripe_df,
      ggplot2::aes(
        xmin = xmin,
        xmax = xmax,
        ymin = -Inf,
        ymax = Inf,
        fill = fill
      ),
      alpha = 0.5,
      inherit.aes = FALSE
    ) +
    ggplot2::scale_fill_identity() +
    ggplot2::geom_point(
      position = ggplot2::position_dodge(width = 0.5)
    ) +
    ggplot2::geom_errorbar(
      ggplot2::aes(
        ymin = lower_ci,
        ymax = upper_ci
      ),
      width = 0.2,
      position = ggplot2::position_dodge(width = 0.5)
    ) +
    ggplot2::scale_color_manual(
      values = c(
        ratio = "#E41A1C",
        numerator = "#377EB8",
        denominator = "#4DAF4A"
      )
    ) +
    ggplot2::scale_x_discrete(drop = FALSE) +
    ggplot2::scale_y_continuous(
      breaks = c(-1, 0, 1, 10, 100),
      labels = function(x) 10^x,
      expand = ggplot2::expansion(mult = c(0.1, 0.1))
    ) +
    ggplot2::geom_hline(
      yintercept = 0.0,
      linetype = "twodash",
      color = "coral2",
      linewidth = 0.5
    ) +
    ggplot2::geom_text(
      ggplot2::aes(
        label = significance,
        y = upper_ci
      ),
      vjust = 0.75,
      hjust = -0.3,
      fontface = "bold",
      position = ggplot2::position_dodge(width = 0.5),
      show.legend = FALSE
    ) +
    ggplot2::labs(
      x = "Cohort",
      y = "Odds Ratio",
      color = "Protein"
    ) +
    ggplot2::facet_wrap(~facet_label, nrow = 1, scales = "free_x") +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1.0),
      panel.spacing.x = grid::unit(1, "lines")
    ) +
    ggplot2::coord_flip(clip = "off")

  return(p)
}




p <- c(ratio_data, individual_data) |>
  dplyr::bind_rows(.id = "cohort") |>
  dplyr::mutate(
    y = factor(y, levels = cohort_order),
    # pair_id = forcats::fct_relevel(pair_id, pair_id_order)
  ) |>
  # Adapt data to the expected inputs of protein_forest_plot()
  # Only retain datapoints for which we have the trio of data for: protein pair, numerator protein and denominator protein
  dplyr::group_by(y, facet_label) |>
  dplyr::filter(dplyr::n_distinct(type) == 3) |>
  dplyr::ungroup() |>
  dplyr::mutate(
    cohort,
    y,
    facet_label,
    x = log10(x),
    lower_ci = log10(ymin),
    upper_ci = log10(ymax),
    colour = type,
    significance = pval_to_stars(p_value),
    .keep = "none"
  ) |>
  .protein_forest_plot()

ggplot2::ggsave(
  plot = p,
  filename = output_figure_file,
  width = 14,
  height = 4,
  units = "in"
)
