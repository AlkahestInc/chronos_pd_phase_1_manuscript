## SPDX-FileCopyrightText: 2026 Alkahest, Inc. <chronospd_data_access@alkahest.com>
##
## SPDX-License-Identifier: MIT

here::i_am("src/40_manuscript/main_figures/figure_3/figure_3g.R")

input_or_chronic_file <- glue::glue("{get_config()$input_dir}/name=OR_phase1_chronic/version=Jan2025/dataset.parquet")
input_jm_comorbidity_file <- glue::glue("{get_config()$input_dir}/name=JM_comorbidity_results/version=Dec2024/dataset.parquet")

output_figure_file <- glue::glue("{get_config()$output_dir}/main_figures/figure_3g.pdf")
dir.create(dirname(output_figure_file), showWarnings = FALSE, recursive = TRUE)


comorbidities_of_interest <- input_or_chronic_file |>
  arrow::open_dataset() |>
  dplyr::filter(
    !diagnosis_code_3 %in% c("G20", "G21"),
    adjp < 0.05
  ) |>
  dplyr::distinct(
    ICD10_Category = diagnosis_code_3
  ) |>
  dplyr::collect()

jm_comorbidity <- input_jm_comorbidity_file |>
  arrow::open_dataset() |>
  dplyr::filter(
    rhat < 1.05,
    p_value < 0.01
  ) |>
  # Restrict to ICD-10 Categories of interest
  dplyr::semi_join(
    comorbidities_of_interest,
    by = dplyr::join_by("diagnosis_code_3" == "ICD10_Category")
  ) |>
  dplyr::mutate(
    # col_label = paste(protein, gene_symbol, sep = ": "),
    # row_label = paste(diagnosis_code_3, ICD10_description, sep = ": "),
    col_label = paste(gene_symbol, protein, sep = "_"),
    row_label = ICD10_description,
    p_value = -log10(dplyr::if_else(
      p_value == 0,
      1e-5,
      p_value
    )),
    .keep = "none"
  ) |>
  dplyr::collect()

logger::log_info("Found {nrow(JM_comorbidity_results_updated)} significant biomarker-comorbidity associations")
logger::log_info("Unique biomarkers: {length(unique(JM_comorbidity_results_updated$protein))}")

# Create a matrix ready for ComplexHeatmap
mat <- jm_comorbidity |>
  tidyr::pivot_wider(
    names_from = col_label,
    values_from = p_value,
    values_fill = 0
  ) |>
  tibble::column_to_rownames("row_label") |>
  as.matrix()

color_fun <- circlize::colorRamp2(c(0, 2, 5), c("white", "yellow", "red"))

hmap <- mat |>
  ComplexHeatmap::Heatmap(
    col = color_fun,
    cluster_rows = TRUE,
    cluster_columns = TRUE,
    row_names_gp = grid::gpar(fontsize = 6),
    column_names_gp = grid::gpar(fontsize = 6),
    show_heatmap_legend = FALSE,

  )

# Create a manual legend object
lgd = ComplexHeatmap::Legend(
  title = "-log10(p-value)",
  direction = "horizontal",
  title_position = "topcenter",
  legend_width = grid::unit(0.15, "npc"),
  col_fun = color_fun,
  at = 0:6,
  title_gp = grid::gpar(fontsize = 8),
  labels_gp = grid::gpar(fontsize = 6)
)

pdf(file = output_figure_file, width = 5, height = 4)
# Draw the heatmap
ComplexHeatmap::draw(
  hmap
)
# Draw the legend manually in the bottom-right corner using viewport coordinates
# x and y are in 'npc' (normalized parent coordinates) from 0 to 1
ComplexHeatmap::draw(
  lgd,
  x = grid::unit(0.80, "npc"),
  y = grid::unit(0.13, "npc"),
  just = c("right", "bottom")
)
dev.off()

logger::log_info("Saved figure to: {output_figure_file}")
