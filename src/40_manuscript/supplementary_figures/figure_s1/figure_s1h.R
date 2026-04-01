## SPDX-FileCopyrightText: 2026 Alkahest, Inc. <chronospd_data_access@alkahest.com>
##
## SPDX-License-Identifier: MIT

here::i_am("src/40_manuscript/supplementary_figures/figure_s1/figure_s1h.R")

chronos_donor_sample_summary_file <- glue::glue("{get_config()$input_dir}/name=chronos_donor_sample_summary/version=v2/dataset.parquet")

output_figure_file <- glue::glue("{get_config()$output_dir}/supplementary_figures/figure_s1/figure_s1h.pdf")
dir.create(dirname(output_figure_file), showWarnings = FALSE, recursive = TRUE)

donor_sample_summary_hashed <- chronos_donor_sample_summary_file |>
  arrow::open_dataset()

plotting_data <- donor_sample_summary_hashed |>
  dplyr::filter(n_sample > 1) |>
  dplyr::select(
    TOKEN_1_hash,
    TOKEN_2_hash,
    donated_for_days,
    n_sample
  ) |>
  dplyr::mutate(
    mean_interval = donated_for_days / (n_sample - 1),
  ) |>
  dplyr::collect()

p <- plotting_data |>
  ggplot2::ggplot(
    ggplot2::aes(x = mean_interval),
  ) +
  ggplot2::geom_histogram(
    binwidth = 1,
    boundary = 0,
  ) +
  ggplot2::scale_y_continuous(
    breaks = seq(0, 400000, by = 50000),
    labels = scales::label_number(scale_cut = scales::cut_short_scale())
  ) +
  ggplot2::scale_x_continuous(
    limits = c(0, 50),
    breaks = seq(0, 50, by = 5),
    labels = scales::number_format(accuracy = 1),
  ) +
  ggplot2::labs(
    title = "Mean Sample Interval per Donor",
    x = "Mean Sample Interval\n(Days)",
    y = "Number of Donors",
  ) +
  ggplot2::annotate(
    "label",
    x = Inf,
    y = Inf,
    label = summary_text(plotting_data, mean_interval),
    family = "mono",
    hjust = 1,
    vjust = 1,
    size = 3.5,
    linewidth = 0,
    label.padding = grid::unit(0.25, "lines"),
    fill = "white",
    alpha = 0.8,
    lineheight = 0.9
  ) +
  ggplot2::theme_minimal() +
  ggplot2::theme(
    panel.grid.minor = ggplot2::element_blank(),
    panel.grid.major = ggplot2::element_blank(),
    axis.ticks = ggplot2::element_line(),
    axis.line = ggplot2::element_line()
  )

ggplot2::ggsave(
  filename = output_figure_file,
  plot = p,
  width = 4.5, height = 3.0, units = "in",
  dpi = 600,
  bg = "white",
)
