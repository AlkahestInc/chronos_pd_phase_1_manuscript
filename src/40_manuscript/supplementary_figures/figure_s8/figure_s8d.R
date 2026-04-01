## SPDX-FileCopyrightText: 2026 Alkahest, Inc. <chronospd_data_access@alkahest.com>
##
## SPDX-License-Identifier: MIT

here::i_am("src/40_manuscript/supplementary_figures/figure_s8/figure_s8d.R")

output_figure_file <- glue::glue("{get_config()$output_dir}/supplementary_figures/figure_s8/figure_s8d.pdf")
dir.create(dirname(output_figure_file), showWarnings = FALSE, recursive = TRUE)

model_group_results_files <- tibble::tribble(
  ~platform,                ~method,           ~file_path,                                                                                                     ~id_col,        ~q_val_col,   ~df_col,  ~coef_col,
  "SomaScan v5",            "1 sample",        glue::glue("{get_config()$input_dir}/name=somalogic_1last_model_group_results/version=none/dataset.parquet"),   "protein_id",   "adj",        "df",     "coef.names",
  "SomaScan v5",            "1-2 samples",     glue::glue("{get_config()$input_dir}/name=somalogic_1to2_model_group_results/version=none/dataset.parquet"),    "ID",           "adj.p",      "DF",     "coef.names",
  "SomaScan v5",            "1-3 samples",     glue::glue("{get_config()$input_dir}/name=somalogic_1to3_model_group_results/version=none/dataset.parquet"),    "ID",           "adj.p",      "DF",     "coef.names",
  "SomaScan v5",            "1-4 samples",     glue::glue("{get_config()$input_dir}/name=somalogic_1to4_model_group_results/version=none/dataset.parquet"),    "ID",           "adj.p",      "DF",     "coef.names",
  "SomaScan v5",            "1-5 samples",     glue::glue("{get_config()$input_dir}/name=somalogic_1to5_model_group_results/version=none/dataset.parquet"),    "ID",           "adj.p",      "DF",     "coef.names",

  "Olink Explore HT",       "1 sample",        glue::glue("{get_config()$input_dir}/name=olink_1last_model_group_results/version=none/dataset.parquet"),       "protein_id",   "adj",        "df",     "coef.names",
  "Olink Explore HT",       "1-2 samples",     glue::glue("{get_config()$input_dir}/name=olink_1to2_model_group_results/version=none/dataset.parquet"),        "ID",           "adj.p",      "DF",     "coef.names",
  "Olink Explore HT",       "1-3 samples",     glue::glue("{get_config()$input_dir}/name=olink_1to3_model_group_results/version=none/dataset.parquet"),        "ID",           "adj.p",      "DF",     "coef.names",
  "Olink Explore HT",       "1-4 samples",     glue::glue("{get_config()$input_dir}/name=olink_1to4_model_group_results/version=none/dataset.parquet"),        "ID",           "adj.p",      "DF",     "coef.names",
  "Olink Explore HT",       "1-5 samples",     glue::glue("{get_config()$input_dir}/name=olink_1to5_model_group_results/version=none/dataset.parquet"),        "ID",           "adj.p",      "DF",     "coef.names",

  "Alamar NULISA CS 120",   "1 sample",        glue::glue("{get_config()$input_dir}/name=alamar_1last_model_group_results/version=none/dataset.parquet"),      "protein_id",   "adj",        "df",     "coef.names",
  "Alamar NULISA CS 120",   "1-2 samples",     glue::glue("{get_config()$input_dir}/name=alamar_1to2_model_group_results/version=none/dataset.parquet"),       "ID",           "adj.p",      "DF",     "coef.names",
  "Alamar NULISA CS 120",   "1-3 samples",     glue::glue("{get_config()$input_dir}/name=alamar_1to3_model_group_results/version=none/dataset.parquet"),       "ID",           "adj.p",      "DF",     "coef.names",
  "Alamar NULISA CS 120",   "1-4 samples",     glue::glue("{get_config()$input_dir}/name=alamar_1to4_model_group_results/version=none/dataset.parquet"),       "ID",           "adj.p",      "DF",     "coef.names",
  "Alamar NULISA CS 120",   "1-5 samples",     glue::glue("{get_config()$input_dir}/name=alamar_1to5_model_group_results/version=none/dataset.parquet"),       "ID",           "adj.p",      "DF",     "coef.names"
)

.read_and_harmonize <- function(file_path, id_col, q_val_col, df_col, coef_col, platform, method) {
  data <- arrow::open_dataset(file_path)

  data <- data |>
    dplyr::rename(
      id = !!rlang::sym(id_col),
      p.adj = !!rlang::sym(q_val_col),
      df = !!rlang::sym(df_col),
      coef_names = !!rlang::sym(coef_col)
    ) |>
    dplyr::mutate(
      platform = platform,
      method = method,
      .before = 1
    ) |>
    dplyr::collect()

  return(data)
}

combined_results_test <- model_group_results_files |>
  purrr::pmap(.read_and_harmonize) |>
  purrr::list_rbind()

plotting_data <- combined_results_test |>
  dplyr::select(
    platform,
    id,
    method,
    p.value,
    p.adj
  ) |>
  tidyr::pivot_longer(
    cols = c(p.value, p.adj),
    names_to = "metric",
    values_to = "value"
  ) |>
  # For each value, determine which thresholds it meets and create multiple rows
  dplyr::mutate(
    thresholds_met = dplyr::case_when(
      value < 0.001 ~ list(c("< 0.05", "< 0.01", "< 0.001")),
      value < 0.01 ~ list(c("< 0.05", "< 0.01")),
      value < 0.05 ~ list(c("< 0.05")),
      .default = list(c("≥ 0.05"))
    )
  ) |>
  tidyr::unnest_longer(thresholds_met) |>
  dplyr::mutate(
    threshold_met = factor(
      thresholds_met,
      levels = c("≥ 0.05", "< 0.05", "< 0.01", "< 0.001")
    ),
    metric = factor(metric, levels = c("p.value", "p.adj"))
  ) |>
  dplyr::select(-thresholds_met)

p <- plotting_data |>
  dplyr::group_by(
    platform,
    method,
    metric,
    threshold_met
  ) |>
  dplyr::summarise(
    n_detect = dplyr::n_distinct(id),
    .groups = "drop"
  ) |>
  dplyr::filter(
    threshold_met != "≥ 0.05"
  ) |>
  ggplot2::ggplot(
    ggplot2::aes(
      x = method,
      y = n_detect,
      group = platform,
      color = platform
    )
  ) +
  ggplot2::geom_line() +
  ggplot2::geom_point() +
  ggplot2::scale_color_manual(
    values = c(
      "Olink Explore HT" = "#00B0F080",
      "Alamar NULISA CS 120" = "#FFB30080",
      "SomaScan v5" = "#00796B80"
    )
  ) +
  ggplot2::scale_y_continuous(
    trans = "log2"
  ) +
  ggplot2::labs(
    x = "",
    y = "Number of significant proteins",
    color = "Platform"
  ) +
  ggplot2::theme_minimal() +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
    legend.title = ggplot2::element_blank(),
    legend.position = "right",
    panel.spacing.y = grid::unit(2, "lines")
  ) +
  ggh4x::facet_grid2(
    metric ~ threshold_met,
    scales = "free_y",
    independent = "y",
    labeller = ggplot2::labeller(
      metric = c(
        "p.value" = "p-value",
        "p.adj" = "adj. p-value"
      )
    )
  )

ggplot2::ggsave(
  filename = output_figure_file,
  plot = p,
  width = 8.0,
  height = 3.0,
  units = "in",
  dpi = 600,
  bg = "white"
)
