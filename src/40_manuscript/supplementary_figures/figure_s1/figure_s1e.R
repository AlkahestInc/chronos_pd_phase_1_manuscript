## SPDX-FileCopyrightText: 2026 Alkahest, Inc. <chronospd_data_access@alkahest.com>
##
## SPDX-License-Identifier: MIT

here::i_am("src/40_manuscript/supplementary_figures/figure_s1/figure_s1e.R")

chronos_donor_demographics_file <- glue::glue("{get_config()$input_dir}/name=chronos_donor_demographics/version=v2/dataset.parquet")

output_figure_file <- glue::glue("{get_config()$output_dir}/supplementary_figures/figure_s1/figure_s1e.pdf")
dir.create(dirname(output_figure_file), showWarnings = FALSE, recursive = TRUE)


plotting_data <- chronos_donor_demographics_file |>
  arrow::open_dataset() |>
  # Limit to donors with samples
  dplyr::filter(
    SAMPLE_COUNT >= 1
  ) |>
  dplyr::group_by(RACE_INDICATOR) |>
  dplyr::tally() |>
  dplyr::collect()

p <- plotting_data |>
  dplyr::rename(
    Ethnic_Background = RACE_INDICATOR,
  ) |>
  dplyr::mutate(
    Ethnic_Background = dplyr::recode_values(
      Ethnic_Background,
      "AFRICAN AMERICAN" ~ "Black",
      "CAUCASIAN" ~ "White",
      "ASIAN" ~ "Asian",
      "HISPANIC/LATINO" ~ "Hispanic/Latino",
      "N/A-DECLINES TO ANSWER" ~ "N/A or Decline to Answer",
      "OTHER" ~ "Other"
    )
  ) |>
  ggplot2::ggplot(
    ggplot2::aes(
      y = forcats::fct_reorder(Ethnic_Background, n),
      x = n,
      fill = Ethnic_Background
    )
  ) +
  ggplot2::geom_bar(
    stat = "identity",
    color = "black"
  ) +
  cowplot::theme_cowplot() +
  ggplot2::geom_text(
    ggplot2::aes(
      label = paste0(round(n / sum(n) * 100, 1), "%"),
      hjust = -0.1
    )
  ) +
  ggplot2::scale_fill_manual(
    values = c(
      "Asian" = "#CD0BBC",
      "Black" = "black",
      "Hispanic/Latino" = "#2297E6",
      "N/A or Decline to Answer" = "#F5C710",
      "Other" = "#DF536B",
      "White" = "white"
    )
  ) +
  ggplot2::theme(
    legend.position = "none",
    legend.box = "horizontal",
    legend.box.just = "center",
    legend.key.width = grid::unit(2, "cm")
  ) +
  ggplot2::scale_x_continuous(
    breaks = seq(0, 1500000, by = 250000),
    expand = ggplot2::expansion(mult = c(0.05, 0.175)),
    labels = scales::label_number(scale_cut = scales::cut_short_scale())
  ) +
  ggplot2::guides(
    color = ggplot2::guide_legend(nrow = 1)
  ) +
  ggplot2::labs(
    x = "Number of Individuals",
    y = "Ethnic Background",
    fill = " "
  )

ggplot2::ggsave(
  filename = output_figure_file,
  plot = p,
  width = 7,
  height = 3,
  units = "in",
  dpi = 600,
  bg = "white"
)
