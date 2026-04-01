## SPDX-FileCopyrightText: 2026 Alkahest, Inc. <chronospd_data_access@alkahest.com>
##
## SPDX-License-Identifier: MIT

here::i_am("src/40_manuscript/supplementary_figures/figure_s1/figure_s1c.R")

chronos_donor_sample_summary_file <- glue::glue("{get_config()$input_dir}/name=chronos_donor_sample_summary/version=v2/dataset.parquet")

output_figure_file <- glue::glue("{get_config()$output_dir}/supplementary_figures/figure_s1/figure_s1c.pdf")
dir.create(dirname(output_figure_file), showWarnings = FALSE, recursive = TRUE)

donor_sample_summary_hashed <- chronos_donor_sample_summary_file |>
  arrow::open_dataset()

p <- donor_sample_summary_hashed |>
  dplyr::select(
    TOKEN_1_hash,
    TOKEN_2_hash,
    current_age:age_stopped,
  ) |>
  dplyr::collect() |>
  tidyr::pivot_longer(
    cols = current_age:age_stopped,
    names_to = "measure",
    values_to = "age",
  ) |>
  dplyr::mutate(
    measure = factor(measure, levels = c("age_started", "age_stopped", "current_age")),
  ) |>
  ggplot2::ggplot(
    ggplot2::aes(
      x = age,
      y = measure,
      fill = measure,
    )
  ) +
  ggridges::stat_density_ridges(
    scale = 1.5,
    rel_min_height = 0.001,
    alpha = 0.6,
  ) +
  ggplot2::scale_x_continuous(
    expand = c(0, 0),
    limits = c(0, 100),
  ) +
  ggplot2::scale_y_discrete(
    expand = c(0, 0),
    labels = c(
      "age_started" = "Age at First Sample",
      "age_stopped" = "Age at Last Sample",
      "current_age" = "Current Age"
    )
  ) +
  ggplot2::coord_cartesian(clip = "off") +
  ggridges::theme_ridges(
    grid = FALSE,
    center_axis_labels = TRUE
  ) +
  ggplot2::theme(
    axis.title.y = ggplot2::element_blank()
  ) +
  ggplot2::labs(
    title = "Chronos",
    subtitle = glue::glue("As of {lubridate::today()}"),
    x = "Age (years)",
  ) +
  ggplot2::scale_fill_manual(
    values = c("#348c43", "#6bcdee", "#4d7bac"),
    guide = "none"
  )

ggplot2::ggsave(
  filename = output_figure_file,
  plot = p,
  width = 4.5, height = 3.0, units = "in",
  dpi = 600,
  bg = "white",
)
