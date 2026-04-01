## SPDX-FileCopyrightText: 2026 Alkahest, Inc. <chronospd_data_access@alkahest.com>
##
## SPDX-License-Identifier: MIT

here::i_am("src/40_manuscript/main_figures/figure_3/figure_3f.R")

input_dataset_files <- list(
  alamar = glue::glue("{get_config()$input_dir}/name=chronos_alamar_joint_model-pd_master_summary/version=29821/dataset.parquet"),
  somalogic = glue::glue("{get_config()$input_dir}/name=chronos_somalogic_joint_model-pd_master_summary/version=29819/dataset.parquet"),
  olink = glue::glue("{get_config()$input_dir}/name=chronos_olink_joint_model-pd_master_summary/version=29820/dataset.parquet")
)

output_figure_file <- glue::glue("{get_config()$output_dir}/main_figures/figure_3f.pdf")
dir.create(dirname(output_figure_file), showWarnings = FALSE, recursive = TRUE)


risk_combined <- purrr::map(
  input_dataset_files,
  ~ {
    arrow::open_dataset(.x) |>
      dplyr::select(
        "probe_id",
        "uniprot_id",
        "gene_symbol",
        "estimate",
        "p_val",
        "q_val",
        "rhat",
        "risk",
        "risk_lower",
        "risk_upper"
      ) |>
      dplyr::collect()
  }
) |>
  dplyr::bind_rows(.id = "dataset") |>
  dplyr::group_by(uniprot_id) |>
  dplyr::slice_min(p_val) |>
  dplyr::ungroup() |>
  dplyr::mutate(
    gene_symbol = forcats::fct_reorder(as.factor(gene_symbol), risk)
  )

plotting_data <- risk_combined |>
  dplyr::filter(
    risk > 0.05 ,
    p_val < 0.005 ,
    rhat < 1.05,
    !gene_symbol %in% c("DDC", "SUMF1", "OIP5")
  )

p <- plotting_data |>
  ggplot2::ggplot(
    ggplot2::aes(
      x = risk,
      y = gene_symbol,
      xmin = risk_lower,
      xmax = risk_upper
    )
  ) +
  ggplot2::geom_point(
    col = "darkblue"
  ) +
  ggplot2::geom_errorbar(
    width = 0.2,
    col = "blue"
  ) +
  ggplot2::scale_x_continuous(
    labels = scales::percent_format(),
    breaks = seq(0, 1, by = 0.05),
    limits = c(0, NA)
  ) +
  ggplot2::labs(
    x = "Risk of developing PD per 0.1 unit change",
    y = "Risk biomarker"
  ) +
  ggplot2::theme_minimal()

ggplot2::ggsave(
  filename = output_figure_file,
  plot = p,
  width = 4.5,
  height = 6,
  units = "in",
  dpi = 600,
  bg = "white"
)
