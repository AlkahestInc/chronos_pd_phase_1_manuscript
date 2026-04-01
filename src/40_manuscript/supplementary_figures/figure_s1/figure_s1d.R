## SPDX-FileCopyrightText: 2026 Alkahest, Inc. <chronospd_data_access@alkahest.com>
##
## SPDX-License-Identifier: MIT

here::i_am("src/40_manuscript/supplementary_figures/figure_s1/figure_s1d.R")

chronos_donor_demographics_file <- glue::glue("{get_config()$input_dir}/name=chronos_donor_demographics/version=v2/dataset.parquet")

output_figure_file <- glue::glue("{get_config()$output_dir}/supplementary_figures/figure_s1/figure_s1d.pdf")
dir.create(dirname(output_figure_file), showWarnings = FALSE, recursive = TRUE)

plotting_data <- chronos_donor_demographics_file |>
  arrow::open_dataset() |>
  # Limit to donors with samples
  dplyr::filter(
    SAMPLE_COUNT >= 1
  ) |>
  dplyr::group_by(GENDER) |>
  dplyr::tally() |>
  dplyr::collect()

p <- plotting_data |>
  dplyr::mutate(
    GENDER = factor(dplyr::recode_values(
      GENDER,
      "F" ~ "Female",
      "M" ~ "Male",
      default = "Multi"
    ),
    levels = c("Female", "Multi", "Male"))
  ) |>
  ggplot2::ggplot(
    ggplot2::aes(
      x = "",
      y = n,
      fill = GENDER,
    )
  ) +
  ggplot2::geom_col(width = 1) +
  ggplot2::coord_polar(theta = "y") +
  ggplot2::scale_fill_manual(values = c(
    Female = "#F8766D",
    Male = "#619CFF"
  )) +
  ggplot2::geom_text(
    ggplot2::aes(
      label = glue::glue("{GENDER}\n{pct}%", pct = round(n / sum(n) * 100, 1))
    ),
    position = ggplot2::position_stack(vjust = 0.5),
  ) +
  ggplot2::labs(title = "Gender Distribution") +
  ggplot2::theme_void() +
  ggplot2::theme(legend.position = "none")

ggplot2::ggsave(
  filename = output_figure_file,
  plot = p,
  width = 3.0, height = 3.0, units = "in",
  dpi = 600,
  bg = "white",
)
