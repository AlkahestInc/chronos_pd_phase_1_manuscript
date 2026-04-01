## SPDX-FileCopyrightText: 2026 Alkahest, Inc. <chronospd_data_access@alkahest.com>
##
## SPDX-License-Identifier: MIT

here::i_am("src/40_manuscript/supplementary_figures/figure_s3/figure_s3a.R")

input_code_category_presence_pre_onset <- glue::glue("{get_config()$input_dir}/name=presence_absence_diagnosis_code_categories_pre_onset/version=30007/dataset.parquet")
input_code_category_presence_post_onset <- glue::glue("{get_config()$input_dir}/name=presence_absence_diagnosis_code_categories_post_onset/version=30008/dataset.parquet")
output_figure_file <- glue::glue("{get_config()$output_dir}/supplementary_figures/figure_s3/figure_s3a.pdf")

# Create output directory if it doesn't exist
dir.create(dirname(output_figure_file), showWarnings = FALSE, recursive = TRUE)

# Load pre-onset data and process into presence-absence matrix for pre-onset
presence_absence_pre_onset <- input_code_category_presence_pre_onset |>
  arrow::read_parquet() |>
  dplyr::mutate(
    value = 1L
  ) |>
  dplyr::distinct(
    PATIENT_ID,
    code_label,
    .keep_all = TRUE
  ) |>
  tidyr::pivot_wider(
    id_cols     = PATIENT_ID,
    names_from  = code_label,
    values_from = value,
    values_fill = list(value = 0)
  )

# Transpose so that patients are columns and codes are rows
presence_absence_pre_onset_t <- presence_absence_pre_onset |>
  tibble::column_to_rownames("PATIENT_ID") |>
  t()

# Load pre-onset data and process into presence-absence matrix for pre-onset
presence_absence_post_onset <- input_code_category_presence_post_onset |>
  arrow::read_parquet() |>
  dplyr::mutate(
    value = 1L
  ) |>
  dplyr::distinct(
    PATIENT_ID,
    code_label,
    .keep_all = TRUE
  ) |>
  tidyr::pivot_wider(
    id_cols     = PATIENT_ID,
    names_from  = code_label,
    values_from = value,
    values_fill = list(value = 0)
  )

# Transpose so that patients are columns and codes are rows
presence_absence_post_onset_t <- presence_absence_post_onset |>
  tibble::column_to_rownames("PATIENT_ID") |>
  t()

# ------------------------------------------------------------------------------
# Create heatmap
# ------------------------------------------------------------------------------

# Get the pre-onset IDs to be in the same order as post-onset IDs
presence_absence_pre_onset_t <- presence_absence_pre_onset_t[, colnames(presence_absence_post_onset_t)]

# Rename patients to distinguish pre-onset and post-onset data
colnames(presence_absence_pre_onset_t) <- paste0(colnames(presence_absence_pre_onset_t), "_pre")
colnames(presence_absence_post_onset_t) <- paste0(colnames(presence_absence_post_onset_t), "_post")

# Filter out "Z" codes for visualization purposes
presence_absence_pre_onset_t <- presence_absence_pre_onset_t[!startsWith(rownames(presence_absence_pre_onset_t), "Z"), ]
presence_absence_post_onset_t <- presence_absence_post_onset_t[!startsWith(rownames(presence_absence_post_onset_t), "Z"), ]

# Filter out NA row
presence_absence_pre_onset_t <- presence_absence_pre_onset_t[
  row.names(presence_absence_pre_onset_t) != "NA", ,
  drop = FALSE
]
presence_absence_post_onset_t <- presence_absence_post_onset_t[
  row.names(presence_absence_post_onset_t) != "NA", ,
  drop = FALSE
]

combined <- cbind(presence_absence_pre_onset_t, presence_absence_post_onset_t)

# Get the column (individual) ordering for plotting - based on post-onset data
post_clust <- combined[, 349:696] |>
  pheatmap::pheatmap(
    cluster_rows = FALSE,
    cluster_cols = TRUE,
    silent       = TRUE
  )
post_clust <- post_clust[["tree_col"]][["order"]]

pre_onset_ind_idx <- post_clust
post_onset_ind_idx <- post_clust + 348

# Save the heatmap to file (silent plotting)
p <- pheatmap::pheatmap(
  combined[, c(pre_onset_ind_idx, post_onset_ind_idx)],
  cellwidth = 1, cellheight = 8,
  border_color = "white",
  color = c("#dddddd", "blue"),
  treeheight_row = 0, treeheight_col = 0,
  show_colnames = FALSE,
  legend_breaks = c(0, 1),
  legend_labels = c("Absent", "Present"),
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  gaps_col = ncol(presence_absence_pre_onset_t),
  filename = NA,
  silent = TRUE
)

ggplot2::ggsave(
  filename = output_figure_file,
  plot = p,
  width = 20, height = 9, units = "in",
  dpi = 600,
  bg = "white"
)
