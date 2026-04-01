## SPDX-FileCopyrightText: 2026 Alkahest, Inc. <chronospd_data_access@alkahest.com>
##
## SPDX-License-Identifier: MIT

here::i_am("src/40_manuscript/main_figures/figure_4/figure_4d.R")

input_comprehensive_data <- glue::glue("{get_config()$input_dir}/name=comprehensive_protein_ratios_data/version=20260128/dataset.parquet")
output_figure_file <- glue::glue("{get_config()$output_dir}/main_figures/figure_4d.pdf")
dir.create(dirname(output_figure_file), showWarnings = FALSE, recursive = TRUE)

# Prepare forest data
forest_data <- input_comprehensive_data |>
  arrow::read_parquet() |>
  # Keep forest plot data: either has odds_ratio values OR is a forest cohort with estimate=NA
  dplyr::filter(!is.na(odds_ratio) | (is.na(estimate) & !is.na(study))) |>
  dplyr::rename(
    protein_pair = ratio,
    cohort = measure,
    cohort_group = study
  ) |>
  dplyr::select(
    protein_pair, gene_symbol1, gene_symbol2, cohort, cohort_group,
    window_info, odds_ratio, or_lower_ci, or_upper_ci, p_value
  )

# ============================================================================
# Create Combined Plot for Selected Ratios
# ============================================================================

# Define the target gene symbol pairs (numerator/denominator)
target_ratios <- tibble::tibble(
  ratio = c(
    "seq.3805.16_o_seq.9191.8",
    "seq.10746.24_o_seq.9191.8",
    "seq.5694.57_o_seq.9191.8",
    "OID45275_o_seq.20579.50",
    "OID45195_o_OID45275",
    "OID45216_o_OID45275",
    "seq.20093.9_o_seq.9191.8",
    "OID43875_o_seq.8275.31",
    "seq.31524.5_o_seq.9191.8",
    "OID45275_o_seq.20579.50",
    "OID44508_o_seq.9191.8"
  ),
)

# Find the protein_pair IDs for these ratios
selected_pairs <- forest_data |>
  dplyr::inner_join(
    target_ratios,
    dplyr::join_by("protein_pair" == "ratio")
  ) |>
  dplyr::pull(protein_pair) |>
  unique()

message(sprintf("\nCreating combined plot for %d selected ratios...\n", length(selected_pairs)))

# Create forest plots for each selected pair
forests <- list()
for (i in seq_along(selected_pairs)) {
  pair <- selected_pairs[i]
  result <- create_forest_plot_for_protein_pair(pair, forest_data)

  if (!is.null(result)) {
    forests[[i]] <- result$plot
  }
}

# Remove NULL entries if any
forests <- Filter(Negate(is.null), forests)

# Define theme modifications
no.legend <- ggplot2::theme(legend.position = "none")
no.yaxis <- ggplot2::theme(
  axis.text.y = ggplot2::element_blank(),
  axis.title.y = ggplot2::element_blank(),
  axis.ticks.y = ggplot2::element_blank()
)

# Apply theme modifications: keep only first y-axis and last legend
for (i in seq_along(forests)) {
  if (i < length(forests)) {
    # Remove legend from all except last plot
    forests[[i]] <- forests[[i]] + no.legend
  }
  if (i > 1) {
    # Remove y-axis from all except first plot
    forests[[i]] <- forests[[i]] + no.yaxis
  }
}

# Combine plots using patchwork in one row
combined_plot <- patchwork::wrap_plots(forests, nrow = 1)

# Save combined plot
pdf(
  output_figure_file,
  width = 16,
  height = 5
)
print(combined_plot)
dev.off()
