## SPDX-FileCopyrightText: 2026 Alkahest, Inc. <chronospd_data_access@alkahest.com>
##
## SPDX-License-Identifier: MIT

here::i_am("src/40_manuscript/main_figures/figure_1/figure_1b.R")

input_counts_dataset_file <- glue::glue("{get_config()$input_dir}/name=disease_counts/version=none/dataset.tsv.gz")

output_figure_file <- glue::glue("{get_config()$output_dir}/main_figures/figure_1b.pdf")
dir.create(dirname(output_figure_file), showWarnings = FALSE, recursive = TRUE)

#####
# Load the count data
#####
counts <- input_counts_dataset_file |>
  readr::read_tsv(
    col_types = readr::cols_only(
      group = readr::col_character(),
      label = readr::col_character(),
      id = readr::col_character(),
      n_donors = readr::col_integer()
    )
  ) |>
  dplyr::mutate(
    label = stringr::str_replace_all(label, stringr::fixed("\\n"), "\n"),
    group = stringr::str_replace_all(group, stringr::fixed("\\n"), "\n"),
    group = factor(group, levels = unique(group))
  )

#####
# Define colour theme for each group
#####
group_theme <- tibble::tribble(
  ~group,                                        ~background_colour,   ~title,        ~datapoint,
  "Oncology",                                    "#E4FAF2",            "#55B98D",     "#55B98D",
  "Neurological\nand\nocular\ndisorders",        "#E9F0F1",            "#5D8890",     "#5D8890",
  "Immunological\nand\nInflammatory\nDiseases",  "#EEF3F5",            "#82B1CB",     "#82B1CB",
  "Cardiovascular\nDiseases",                    "#EFF3F9",            "#2F65D8",     "#2F65D8",
  "Infectious\nDiseases",                        "#F5F9F9",            "#7FC7C6",     "#7FC7C6",
  "Respiratory\nDiseases",                       "#E2E2E2",            "#171717",     "#171717",
  "Metabolic\nDisorders",                        "#EBFBFC",            "#194245",     "#194245",
)

#####
# Generate Plot
#####
p <- counts |>
  plot_disease_map(
    group_theme = group_theme,
    y = "n_donors",
    y_tick_values = 10L^c(2L:6L)
  )

width <- ((counts |> nrow() * 1.0) + (nlevels(counts$group) * 2 * 1.5)) * 0.2
ggplot2::ggsave(
  filename = output_figure_file,
  plot = p,
  width = width, height = width / 2, units = "in",
  dpi = 600,
  bg = "white"
)
