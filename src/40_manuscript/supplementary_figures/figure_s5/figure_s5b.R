## SPDX-FileCopyrightText: 2026 Alkahest, Inc. <chronospd_data_access@alkahest.com>
##
## SPDX-License-Identifier: MIT

#!/usr/bin/env Rscript
#' ---
#' title: "Supplementary Figure S5b: Biomarker trajectories (3 panels)"
#' description: >
#'   Biomarker trajectories split into three panels showing different temporal groups:
#'   Panel 1 (left): Inflection years before onset (prodromal) - 47 proteoforms
#'   Panel 2 (middle): Early linear change (transient_preclinical) - 12 proteoforms
#'   Panel 3 (right): Late linear and nonlinear change (group3) - 36 proteoforms
#'   Ported from Figure3_PD_biomarker_trajectory_3groups_plot.R (lines 136-186)
#' output: PDF figure (3 panels side by side)
#' ---

# Setup ========================================================================
here::i_am("src/40_manuscript/supplementary_figures/figure_s5/figure_s5b.R")

# Setup logging
logger::log_layout(logger::layout_glue_generator(format = "{level} [{time}] {msg}"))
logger::log_threshold(logger::INFO)

# Check required packages
required_packages <- c(
  "dplyr", "ggplot2", "cowplot", "arrow", "glue", "logger"
)
missing_packages <- required_packages[!sapply(required_packages, requireNamespace, quietly = TRUE)]
if (length(missing_packages) > 0) {
  stop("Missing required packages: ", paste(missing_packages, collapse = ", "))
}

# Configuration ================================================================

# Input data file (same as Figure 3e)
input_data_file <- glue::glue("{get_config()$input_dir}/name=biomarker-combined_predicted-results/version=none/dataset.parquet")

# Output figure file
output_figure_file <- glue::glue("{get_config()$output_dir}/supplementary_figures/figure_s5/figure_s5b.pdf")
dir.create(dirname(output_figure_file), showWarnings = FALSE, recursive = TRUE)

# Define "%nin%" operator for filtering
`%nin%` <- Negate(`%in%`)

# Load combined predicted results
logger::log_info("Loading biomarker trajectory data from: {input_data_file}")
combined_predicted_results <- arrow::open_dataset(input_data_file) |>
  dplyr::collect()

# Re-normalization to handle cross-zero values (same as Figure 3e)
predicted_results_updated <- combined_predicted_results |>
  dplyr::group_by(probe_id) |>
  dplyr::mutate(
    pred_min = min(normalized),
    pred_max = max(normalized),
    zero_normalized = dplyr::case_when(
      sign(pred_min) != sign(pred_max) & abs(pred_min) > abs(pred_max) ~ normalized - abs(pred_max),
      sign(pred_min) != sign(pred_max) & abs(pred_min) < abs(pred_max) ~ normalized + abs(pred_min),
      .default = normalized
    ),
    biomarker_type = factor(
      biomarker_type,
      levels = c("transient_preclinical", "prodromal", "transitional", "progression", "postDx")
    )
  ) |>
  dplyr::filter(
    # Biomarkers affected by medications (exclude from visualization)
    !gene_symbol %in% c("DDC", "NAA10", "OIP5", "ZC3H12C", "CYLD", "TUBA1A")
  )

# Classify biomarker patterns (same as Figure 3e)
biomarker_groups_combined <- predicted_results_updated |>
  dplyr::group_by(gene_symbol) |>
  dplyr::mutate(
    pattern = dplyr::case_when(
      # Group 1: Transient preclinical
      biomarker_type == "transient_preclinical" ~
        dplyr::if_else(
          probe_id %in% c("seq.2515.14"),
          "increased",
          "decreased"
        ),
      # Group 2: Prodromal
      biomarker_type == "prodromal" ~
        dplyr::if_else(
          min(zero_normalized) >= 0,
          "increased",
          "decreased"
        ),
      # Group 3: Late biomarkers (transitional)
      biomarker_type == "transitional" ~
        dplyr::if_else(
          gene_symbol %in% c("TUBA1A", "SUMF1"),
          "fast_increased",
          "fast_decreased"
        ),
      # Group 3: Late biomarkers (Progression)
      biomarker_type == "progression" ~
        dplyr::if_else(
          gene_symbol %in% c("MTAP", "DDC"),
          "fast_increased",
          "fast_decreased"
        ),
      # Group 3: Late biomarkers (PostDx)
      biomarker_type == "postDx" ~
        dplyr::if_else(
          min(zero_normalized) >= 0,
          "slow_increased",
          "slow_decreased"
        ),
      .default = NA_character_
    )
  ) |>
  dplyr::ungroup() |>
  dplyr::mutate(
    normalized_prediction = dplyr::case_when(
      biomarker_type == "transient_preclinical" ~ normalized,
      biomarker_type == "prodromal" ~ zero_normalized,
      biomarker_type == "transitional" ~ zero_normalized,
      biomarker_type == "progression" ~ zero_normalized,
      biomarker_type == "postDx" ~ zero_normalized,
      .default = NA_real_
    ),
    group = dplyr::case_when(
      biomarker_type == "transient_preclinical" ~ "group1",
      biomarker_type == "prodromal" ~ "group2",
      biomarker_type == "transitional" ~ "group3",
      biomarker_type == "progression" ~ "group3",
      biomarker_type == "postDx" ~ "group3",
      .default = NA_character_
    ),
    pattern_groups = paste(group, pattern, sep = "_")
  ) |>
  dplyr::select(
    -normalized,
    -zero_normalized
  )

# Count proteoforms per group
group1_count <- biomarker_groups_combined |>
  dplyr::filter(group == "group1") |>
  dplyr::distinct(probe_id) |>
  nrow()

group2_count <- biomarker_groups_combined |>
  dplyr::filter(group == "group2") |>
  dplyr::distinct(probe_id) |>
  nrow()

group3_count <- biomarker_groups_combined |>
  dplyr::filter(group == "group3") |>
  dplyr::distinct(probe_id) |>
  nrow()

logger::log_info("Group 1 (inflection): {group1_count} proteoforms")
logger::log_info("Group 2 (early): {group2_count} proteoforms")
logger::log_info("Group 3 (late): {group3_count} proteoforms")

# Create three separate panels ================================================
p2 <- biomarker_groups_combined |>
  dplyr::filter(group == "group1") |>
  ggplot2::ggplot(
    ggplot2::aes(
      x = years_since_pd_onset,
      y = normalized_prediction,
      group = probe_id,
      col = pattern_groups
    )
  ) +
  ggplot2::geom_line(alpha = 0.3, lty = 1) +
  ggplot2::geom_smooth(
    ggplot2::aes(
      group = pattern_groups
    ),
    linewidth = 1.5,
    method = "gam",
    formula = y ~ s(x, bs = "cs")
  ) +
  ggplot2::geom_hline(yintercept = 0, lty = 2, col = "black") +
  ggplot2::scale_color_manual(
    name = "pattern",
    values = c(
      "group1_increased" = "blueviolet",
      "group1_decreased" = "violet"
    ),
    labels = c(
      "group1_increased" = "increased",
      "group1_decreased" = "decreased"
    )
  ) +
  ggplot2::guides(color = ggplot2::guide_legend(ncol = 2)) +
  ggplot2::labs(
    title = glue::glue("Inflection years before onset\n{group1_count} proteoforms"),
    y = "Control Normalized Mean Prediction",
    x = "Years Since PD Onset"
  ) +
  ggplot2::theme_minimal() +
  ggplot2::theme(
    panel.grid.minor = ggplot2::element_blank(),
    legend.position = "bottom",
    plot.title = ggplot2::element_text(hjust = 0.5)
  )

p1 <- biomarker_groups_combined |>
  dplyr::filter(group == "group2") |>
  ggplot2::ggplot(
    ggplot2::aes(
      x = years_since_pd_onset,
      y = normalized_prediction,
      group = probe_id,
      col = pattern_groups
    )
  ) +
  ggplot2::geom_line(alpha = 0.3, lty = 1) +
  ggplot2::geom_smooth(
    ggplot2::aes(
      group = pattern_groups
    ),
    linewidth = 1.5,
    method = "gam",
    formula = y ~ s(x, bs = "cs")
  ) +
  ggplot2::geom_hline(yintercept = 0, lty = 2, col = "black") +
  ggplot2::scale_color_manual(
    name = "pattern",
    values = c(
      "group2_increased" = "royalblue",
      "group2_decreased" = "cyan"
    ),
    labels = c(
      "group2_increased" = "increased",
      "group2_decreased" = "decreased"
    )
  ) +
  ggplot2::guides(color = ggplot2::guide_legend(ncol = 2)) +
  ggplot2::labs(
    title = glue::glue("Early linear change\n(10+ years before onset)\n{group2_count} proteoforms"),
    y = "Control Normalized Mean Prediction",
    x = "Years Since PD Onset"
  ) +
  ggplot2::theme_minimal() +
  ggplot2::theme(
    panel.grid.minor = ggplot2::element_blank(),
    legend.position = "bottom",
    plot.title = ggplot2::element_text(hjust = 0.5)
  )

p3 <- biomarker_groups_combined |>
  dplyr::filter(group == "group3") |>
  ggplot2::ggplot(
    ggplot2::aes(
      x = years_since_pd_onset,
      y = normalized_prediction,
      group = probe_id,
      col = pattern_groups
    )
  ) +
  ggplot2::geom_line(alpha = 0.3, lty = 1) +
  ggplot2::geom_smooth(
    ggplot2::aes(
      group = pattern_groups
    ),
    linewidth = 1.5,
    method = "gam",
    formula = y ~ s(x, bs = "cs")
    ) +
  ggplot2::geom_hline(yintercept = 0, lty = 2, col = "black") +
  ggplot2::scale_color_manual(
    name = "pattern",
    values = c(
      "group3_fast_increased" = "red",
      "group3_fast_decreased" = "hotpink",
      "group3_slow_increased" = "darkred",
      "group3_slow_decreased" = "pink"
    ),
    labels = c(
      "group3_fast_increased" = "fast_increased",
      "group3_fast_decreased" = "fast_decreased",
      "group3_slow_increased" = "slow_increased",
      "group3_slow_decreased" = "slow_decreased"
    )
  ) +
  ggplot2::guides(color = ggplot2::guide_legend(ncol = 2)) +
  ggplot2::labs(
    title = glue::glue("Late linear and nonlinear change\n(mostly peri / post onset)\n{group3_count} proteoforms"),
    y = "Control Normalized Mean Prediction",
    x = "Years Since PD Onset"
  ) +
  ggplot2::theme_minimal() +
  ggplot2::theme(
    panel.grid.minor = ggplot2::element_blank(),
    legend.position = "bottom",
    plot.title = ggplot2::element_text(hjust = 0.5)
  )

# Combine three panels side by side: prodromal (left), transient (middle), group3 (right)
# Equal widths now that all legends use 2-column layout with smaller text
p_combined <- cowplot::plot_grid(
  p1, p2, p3,
  ncol = 3,
  align = "h",
  rel_widths = c(1, 1, 1)
)

# Save plot
ggplot2::ggsave(
  filename = output_figure_file,
  plot = p_combined,
  width = 11,
  height = 5.0,
  units = "in",
  dpi = 600,
  bg = "white"
)

logger::log_info("Saved figure to: {output_figure_file}")
