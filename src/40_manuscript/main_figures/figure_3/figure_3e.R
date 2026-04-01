## SPDX-FileCopyrightText: 2026 Alkahest, Inc. <chronospd_data_access@alkahest.com>
##
## SPDX-License-Identifier: MIT

here::i_am("src/40_manuscript/main_figures/figure_3/figure_3e.R")

# ==============================================================================
# Figure 3e: Biomarker trajectories
# ==============================================================================
# Description: Biomarkers identified in Chronos PD phase 1 dataset based on
#              linear, non-linear mixed models, joint modeling of protein levels
#              and time to PD onset. Trajectories are visualized as
#              control-normalized estimates from 12 years before to 9 years
#              after PD onset, grouped into three biomarker types with 8 patterns.
# ==============================================================================

# Input data file (combined biomarker trajectory predictions from joint modeling)
input_data_file <- glue::glue("{get_config()$input_dir}/name=biomarker-combined_predicted-results/version=none/dataset.parquet")

# Output figure file
output_figure_file <- glue::glue("{get_config()$output_dir}/main_figures/figure_3e.pdf")
dir.create(dirname(output_figure_file), showWarnings = FALSE, recursive = TRUE)

# Define "%nin%" operator for filtering
`%nin%` <- Negate(`%in%`)

# Load combined predicted results
logger::log_info("Loading biomarker trajectory data from: {input_data_file}")
combined_predicted_results <- arrow::open_dataset(input_data_file) |>
  dplyr::collect()

# Re-normalization to handle cross-zero values (for visualization)
# If min and max are not in same direction, follow direction of max values
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


biomarker_groups_combined <- predicted_results_updated |>
  dplyr::group_by(gene_symbol) |>
  dplyr::mutate(
    # Classify biomarker patterns according to the biomarker type
    pattern = dplyr::case_when(
      # Group 1: Transient preclinical
      biomarker_type == "transient_preclinical" ~
        dplyr::if_else(
          probe_id %in% c("seq.2515.14"),
          "increased",
          "decreased"
        ),
      # Group 2:Prodromal (grouped logic)
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
    pattern_groups = paste(group, pattern, sep = "_"),
  ) |>
  dplyr::select(
    -normalized,
    -zero_normalized,
  )


# Create combined plot showing all three biomarker groups
p <- biomarker_groups_combined |>
  dplyr::select(
    years_since_pd_onset,
    normalized_prediction,
    probe_id,
    pattern_groups
  ) |>
  ggplot2::ggplot(
    ggplot2::aes(
      x = years_since_pd_onset,
      y = normalized_prediction,
      group = probe_id,
      col = pattern_groups
    )
  ) +
  ggplot2::geom_line(alpha = 0.8, lty = 1) +
  ggplot2::geom_hline(yintercept = 0, lty = 2, col = "black") +
  ggplot2::scale_color_manual(
    values = c(
      "group1_increased" = "blueviolet",
      "group1_decreased" = "violet",
      "group2_increased" = "royalblue",
      "group2_decreased" = "cyan",
      "group3_fast_increased" = "red",
      "group3_fast_decreased" = "hotpink",
      "group3_slow_increased" = "darkred",
      "group3_slow_decreased" = "pink"
    ),
    guide = "none"
  ) +
  ggplot2::labs(
    y = "Control normalized mean prediction",
    x = "Years since PD onset"
  ) +
  ggplot2::theme_minimal() +
  ggplot2::theme(
    axis.text = ggplot2::element_text(size = 12),
    axis.title = ggplot2::element_text(size = 12, face = "bold"),
    legend.position = "right"
  )

# Save plot
ggplot2::ggsave(
  filename = output_figure_file,
  plot = p,
  width = 5,
  height = 5,
  units = "in",
  dpi = 600,
  bg = "white"
)

logger::log_info("Saved figure to: {output_figure_file}")
