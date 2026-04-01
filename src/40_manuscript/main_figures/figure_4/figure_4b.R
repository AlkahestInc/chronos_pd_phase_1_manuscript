## SPDX-FileCopyrightText: 2026 Alkahest, Inc. <chronospd_data_access@alkahest.com>
##
## SPDX-License-Identifier: MIT

here::i_am("src/40_manuscript/main_figures/figure_4/figure_4b.R")

infile <- glue::glue("{get_config()$input_dir}/name=comprehensive_protein_ratios_data/version=20260128/dataset.parquet")
output_figure_file <- glue::glue("{get_config()$output_dir}/main_figures/figure_4b.pdf")
dir.create(dirname(output_figure_file), recursive = TRUE, showWarnings = FALSE)

comprehensive_data <- infile |>
  arrow::open_dataset() |>
  dplyr::collect()

# ============================================================================
# Reshape Data from Long to Wide Format
# ============================================================================

# Filter for rows with estimate values (heatmap data)
heatmap_long <- comprehensive_data |>
  dplyr::filter(!is.na(estimate))

# Pivot estimates to wide format
estimates_wide <- heatmap_long |>
  dplyr::select(ratio, gene_symbol1, gene_symbol2, measure, estimate) |>
  tidyr::pivot_wider(
    names_from = measure,
    values_from = estimate,
    id_cols = c(ratio, gene_symbol1, gene_symbol2)
  )

# Pivot p-values to wide format
pvalues_wide <- heatmap_long |>
  dplyr::select(ratio, gene_symbol1, gene_symbol2, measure, p_value) |>
  tidyr::pivot_wider(
    names_from = measure,
    values_from = p_value,
    id_cols = c(ratio, gene_symbol1, gene_symbol2)
  )

# 1. Chronos windows (ascending by window_start)
chronos_cols <- names(estimates_wide)[grep("^Chronos_", names(estimates_wide))]
chronos_cols <- chronos_cols[order(as.numeric(sub("Chronos_(-?\\d+)_.*", "\\1", chronos_cols)))]

# 2. UKB windows (ascending by window_start: -9, then -5)
ukb_cols <- names(estimates_wide)[grep("^UKB_", names(estimates_wide))]
ukb_cols <- ukb_cols[order(as.numeric(sub("UKB_(-?\\d+)_.*", "\\1", ukb_cols)))]

# 3. PPMI
ppmi_col <- "PPMI"

# 4. GNPC cohorts (C, J, L)
gnpc_cols <- c("GNPC C", "GNPC J", "GNPC L")
gnpc_cols <- gnpc_cols[gnpc_cols %in% names(estimates_wide)]

column_order <- c("ratio", "gene_symbol1", "gene_symbol2", chronos_cols, ukb_cols, ppmi_col, gnpc_cols)
estimates_wide <- estimates_wide |> dplyr::select(dplyr::all_of(column_order))
pvalues_wide <- pvalues_wide |> dplyr::select(dplyr::all_of(column_order))

# ============================================================================
# Prepare Matrices for Heatmap
# ============================================================================

# Create estimate matrix (excluding identifiers)
estimates_mat <- estimates_wide |>
  dplyr::select(-c(ratio, gene_symbol1, gene_symbol2)) |>
  as.matrix()

# Create p-value matrix with significance stars
pvalues_mat <- pvalues_wide |>
  dplyr::mutate(
    dplyr::across(
      -c(ratio, gene_symbol1, gene_symbol2),
      ~ dplyr::case_when(
        is.na(.x) ~ "",
        .x < 0.001 ~ "***",
        .x < 0.01 ~ "**",
        .x < 0.05 ~ "*",
        TRUE ~ ""
      )
    )
  ) |>
  dplyr::select(-c(ratio, gene_symbol1, gene_symbol2)) |>
  as.matrix()

# Set row names to gene symbols
rownames(estimates_mat) <- paste(estimates_wide$gene_symbol1, estimates_wide$gene_symbol2, sep = " / ")
rownames(pvalues_mat) <- paste(pvalues_wide$gene_symbol1, pvalues_wide$gene_symbol2, sep = " / ")

# ============================================================================
# Define Column Groups for Visual Separation
# ============================================================================

n_cols <- ncol(estimates_mat)

# Define column groups based on data structure
# Chronos: 8 windows, UKB: 2 windows, PPMI: 1, GNPC: remaining
column_groups <- c(
  rep("Group1", 8), # Chronos windows
  rep("Group2", 2), # UKB windows
  rep("Group3", 1), # PPMI
  rep("Group4", n_cols - 11) # GNPC cohorts
)

column_groups <- factor(column_groups,
  levels = c("Group1", "Group2", "Group3", "Group4")
)

# ============================================================================
# STEP 5: Create Heatmap
# ============================================================================

# Define color palette for estimates (log odds)
col_fun <- circlize::colorRamp2(
  c(-2, 0, 2),
  c("blue", "white", "red")
)

# Create heatmap with ComplexHeatmap
heatmap_plot <- ComplexHeatmap::Heatmap(
  estimates_mat,
  name = "Estimate\n(log odds)",
  col = col_fun,

  # Add significance stars
  cell_fun = function(j, i, x, y, width, height, fill) {
    if (!is.na(pvalues_mat[i, j]) && pvalues_mat[i, j] != "") {
      grid::grid.text(
        pvalues_mat[i, j],
        x,
        y - grid::unit(1, "mm"),
        gp = grid::gpar(fontsize = 10, col = "black")
      )
    }
  },

  # Clustering options
  cluster_rows = TRUE,
  cluster_columns = FALSE,

  # Column grouping and separation
  column_split = column_groups,
  column_title = NULL,
  column_gap = grid::unit(3, "mm"),

  # Row and column styling
  show_row_names = TRUE,
  row_names_gp = grid::gpar(fontsize = 8),
  column_names_rot = 45,
  column_names_gp = grid::gpar(fontsize = 9),

  # NA handling
  na_col = "grey90",

  # Legend
  heatmap_legend_param = list(
    title = "Estimate\n(log odds)",
    title_gp = grid::gpar(fontsize = 10, fontface = "bold"),
    labels_gp = grid::gpar(fontsize = 9),
    grid_height = grid::unit(4, "mm"),
    grid_width = grid::unit(4, "mm")
  )
)

# ============================================================================
# Save Heatmap
# ============================================================================
pdf(
  output_figure_file,
  width = 18,
  height = 12
)
ComplexHeatmap::draw(heatmap_plot)
dev.off()
