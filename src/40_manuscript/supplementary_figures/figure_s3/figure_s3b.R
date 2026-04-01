## SPDX-FileCopyrightText: 2026 Alkahest, Inc. <chronospd_data_access@alkahest.com>
##
## SPDX-License-Identifier: MIT

here::i_am("src/40_manuscript/supplementary_figures/figure_s3/figure_s3b.R")

# Define input files
chronos_input_file <- glue::glue("{get_config()$input_dir}/name=chronos-pd_phase_1_odds_ratio_results/version=none/dataset.parquet")
us_input_file <- glue::glue("{get_config()$input_dir}/name=us_population_high_confidence_odds_ratio_results/version=none/dataset.parquet")

# Define output file
output_figure_file <- glue::glue("{get_config()$output_dir}/supplementary_figures/figure_s3b.pdf")
dir.create(dirname(output_figure_file), showWarnings = FALSE, recursive = TRUE)

# Load inputs and separate results into pre- and post-onset
chronos_pre <- chronos_input_file |>
  arrow::open_dataset() |>
  dplyr::filter(
    stage == "before"
  ) |>
  dplyr::select(
    DIAGNOSIS_CODE_CATEGORY,
    PValue,
    OR,
    logOR,
    FDR,
    label
  ) |>
  dplyr::collect()

us_pre <- us_input_file |>
  arrow::open_dataset() |>
  dplyr::filter(
    stage == "before"
  ) |>
  dplyr::select(
    DIAGNOSIS_CODE_CATEGORY,
    PValue,
    OR,
    logOR,
    FDR,
    label
  ) |>
  dplyr::collect()

chronos_post <- chronos_input_file |>
  arrow::open_dataset() |>
  dplyr::filter(
    stage == "after"
  ) |>
  dplyr::select(
    DIAGNOSIS_CODE_CATEGORY,
    PValue,
    OR,
    logOR,
    FDR,
    label
  ) |>
  dplyr::collect()


us_post <- us_input_file |>
  arrow::open_dataset() |>
  dplyr::filter(
    stage == "after"
  ) |>
  dplyr::select(
    DIAGNOSIS_CODE_CATEGORY,
    PValue,
    OR,
    logOR,
    FDR,
    label
  ) |>
  dplyr::collect()

# Identify codes in common between Chronos and US - pre onset
pre_overlapping_codes <- chronos_pre |>
  dplyr::inner_join(
    us_pre,
    by = c("DIAGNOSIS_CODE_CATEGORY"),
    suffix = c("_phase1", "_uspop")
  ) |>
  dplyr::mutate(
    sig_phase1 = FDR_phase1 < 0.05,
    sig_uspop = FDR_uspop < 0.05
  )

# Identify codes in common between Chronos and US - pre onset
post_overlapping_codes <- chronos_post |>
  dplyr::inner_join(
    us_post,
    by = c("DIAGNOSIS_CODE_CATEGORY"),
    suffix = c("_phase1", "_uspop")
  ) |>
  dplyr::mutate(
    sig_phase1 = FDR_phase1 < 0.05,
    sig_uspop = FDR_uspop < 0.05
  )

plotting_data_pre <- pre_overlapping_codes |>
  dplyr::mutate(
    sig_category = dplyr::case_when(
      sig_phase1 & sig_uspop ~ "Significant in both",
      sig_uspop & !sig_phase1 ~ "Significant only in US Population-Level PD",
      sig_phase1 & !sig_uspop ~ "Significant only in Chronos-PD Phase 1",
      TRUE ~ "Not significant"
    ),
    direction_agreement = ifelse(
      (OR_phase1 >= 1 & OR_uspop >= 1) | (OR_phase1 < 1 & OR_uspop < 1),
      "Same direction", "Opposite direction"
    )
  ) |>
  # Filter out extreme values for visualisation purposes
  dplyr::filter(
    logOR_phase1 > -3.5,
    logOR_uspop > -3.5,
    logOR_phase1 < 3.5,
    logOR_uspop < 3.5
  ) |>
  dplyr::mutate(
    sig_category = factor(
      sig_category,
      levels = c(
        "Not significant",
        "Significant only in US Population-Level PD",
        "Significant only in Chronos-PD Phase 1",
        "Significant in both"
      )
    ),
    on_top = sig_category == "Significant in both"
  ) |>
  dplyr::arrange(on_top) # FALSE rows first, TRUE rows last

plotting_data_post <- post_overlapping_codes |>
  dplyr::mutate(
    sig_category = dplyr::case_when(
      sig_phase1 & sig_uspop ~ "Significant in both",
      sig_uspop & !sig_phase1 ~ "Significant only in US Population-Level PD",
      sig_phase1 & !sig_uspop ~ "Significant only in Chronos-PD Phase 1",
      TRUE ~ "Not significant"
    ),
    direction_agreement = ifelse(
      (OR_phase1 >= 1 & OR_uspop >= 1) | (OR_phase1 < 1 & OR_uspop < 1),
      "Same direction", "Opposite direction"
    )
  ) |>
  # Filter out extreme values for visualisation purposes
  dplyr::filter(
    logOR_phase1 > -3.5,
    logOR_uspop > -3.5,
    logOR_phase1 < 3.5,
    logOR_uspop < 3.5
  ) |>
  dplyr::mutate(
    sig_category = factor(
      sig_category,
      levels = c(
        "Not significant",
        "Significant only in US Population-Level PD",
        "Significant only in Chronos-PD Phase 1",
        "Significant in both"
      )
    ),
    on_top = sig_category == "Significant in both"
  ) |>
  dplyr::arrange(on_top) # FALSE rows first, TRUE rows last

plotting_data_combined <- dplyr::bind_rows(
  plotting_data_pre,
  plotting_data_post,
  .id = "Stage"
) |>
  dplyr::mutate(
    Stage = dplyr::case_when(
      Stage == 1 ~ "Pre-Onset",
      Stage == 2 ~ "Post-Onset"
    ),
    Stage = factor(Stage, levels = c("Pre-Onset", "Post-Onset"))
  )

p <- ggplot2::ggplot(
  plotting_data_combined,
  ggplot2::aes(x = OR_phase1, y = OR_uspop)
) +
  ggplot2::geom_point(
    ggplot2::aes(color = sig_category, shape = direction_agreement),
    alpha = 0.7, size = 2
  ) +
  ggplot2::geom_vline(xintercept = 1, linetype = "dashed", color = "blue") +
  ggplot2::geom_hline(yintercept = 1, linetype = "dashed", color = "blue") +
  ggplot2::scale_x_log10() +
  ggplot2::scale_y_log10() +
  ggplot2::scale_color_manual(values = c(
    "Significant in both" = "red",
    "Significant only in US Population-Level PD" = "#009E73",
    "Significant only in Chronos-PD Phase 1" = "#F0E442",
    "Not significant" = "grey70"
  )) +
  ggplot2::scale_shape_manual(
    values = c("Same direction" = 16, "Opposite direction" = 4)
  ) +
  ggplot2::labs(
    x = "Odds Ratio (Chronos-PD Phase 1 - log scale)",
    y = "Odds Ratio (US Population-Level PD - log scale)", 
    color = "Significance Category",
    shape = "Concordance",
    title = "Directional Concordance of Odds Ratios"
  ) +
  ggplot2::theme_minimal(base_size = 12) +
  ggplot2::facet_wrap(~Stage)

ggplot2::ggsave(
  plot = p,
  filename = output_figure_file,
  width = 9,
  height = 4,
  units = "in"
)
