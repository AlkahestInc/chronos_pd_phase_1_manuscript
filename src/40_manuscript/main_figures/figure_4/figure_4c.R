## SPDX-FileCopyrightText: 2026 Alkahest, Inc. <chronospd_data_access@alkahest.com>
##
## SPDX-License-Identifier: MIT

here::i_am("src/40_manuscript/main_figures/figure_4/figure_4c.R")

# Input data file (combined protein ratio results from all cohorts)
input_data_file <- glue::glue("{get_config()$input_dir}/name=protein_ratio-statistics-summary/version=20260128/dataset.rds")

# Output figure file
output_figure_file <- glue::glue("{get_config()$output_dir}/main_figures/figure_4c.pdf")
dir.create(dirname(output_figure_file), showWarnings = FALSE, recursive = TRUE)

# Load combined protein ratio data
logger::log_info("Loading combined protein ratio data from: {input_data_file}")
protein_ratio_data <- readr::read_rds(input_data_file)

# Extract p-values for validation cohorts
pval_data <- protein_ratio_data$p.values

# For UKB, find minimum p-value across time windows
ukb_min_pval <- apply(
  data.frame(pval_data$`UKB_-9_-5`, pval_data$`UKB_-5_-1`),
  1,
  function(x) min(x, na.rm = TRUE)
)

# Replace Inf (when NA for both windows) with NA
ukb_min_pval[is.infinite(ukb_min_pval)] <- NA

logger::log_info("Found {sum(is.infinite(apply(data.frame(pval_data$`UKB_-9_-5`, pval_data$`UKB_-5_-1`), 1, function(x) min(x, na.rm = TRUE))))} protein ratios with NA in both UKB windows")

# Create summary dataframe with validation results across cohorts
ratio_summary <- data.frame(
  gene_symbol1 = pval_data$gene_symbol1,
  gene_symbol2 = pval_data$gene_symbol2,
  GNPC_C = pval_data$`GNPC C_p_value`,
  GNPC_J = pval_data$`GNPC J_p_value`,
  GNPC_L = pval_data$`GNPC L_p_value`,
  PPMI = pval_data$PPMI_p_value,
  UKB_min = ukb_min_pval
)

# Build edge list with weights = number of cohorts where ratio validated (p < 0.05)
edge_list <- data.frame(
  from = ratio_summary$gene_symbol1,
  to = ratio_summary$gene_symbol2,
  weight = apply(
    ratio_summary[, -c(1, 2)],
    1,
    function(x) sum(x < 0.05, na.rm = TRUE)
  )
)

# Sort by weight (number of validations) in descending order
edge_list <- edge_list[order(edge_list$weight, decreasing = TRUE), ]

# Handle duplicate protein pairs (A-B same as B-A)
# Sort protein names alphabetically within each pair
pair_ids <- vapply(
  seq_len(nrow(edge_list)),
  function(i) {
    sorted_pair <- sort(c(edge_list$from[i], edge_list$to[i]))
    paste(sorted_pair[1], sorted_pair[2], sep = "_")
  },
  character(1)
)

logger::log_info("Removing {sum(duplicated(pair_ids))} duplicate protein pairs")

# Remove duplicates (keeping the one with higher weight due to sorting)
edge_list <- edge_list[!duplicated(pair_ids), ]

# Remove rows with NA gene symbols
edge_list <- na.omit(edge_list)

logger::log_info("Total unique protein pairs: {nrow(edge_list)}")

# Select protein ratios validated in at least 2 independent cohorts
edge_list_filtered <- edge_list[edge_list$weight >= 2, ]

logger::log_info("Protein pairs validated in ≥2 cohorts: {nrow(edge_list_filtered)}")

# Create network graph
network_graph <- igraph::graph_from_data_frame(edge_list_filtered, directed = FALSE)

logger::log_info("Network has {igraph::vcount(network_graph)} nodes and {igraph::ecount(network_graph)} edges")

# Calculate node properties
node_degree <- igraph::degree(network_graph)
log_degree <- log(node_degree + 1)
norm_log_degree <- (log_degree - min(log_degree)) / (max(log_degree) - min(log_degree))
contrast_norm <- norm_log_degree^0.3

# Node colors: light to dark blue based on connection degree
node_palette <- grDevices::colorRampPalette(c("#deebf7", "#08306b"))(100)
vertex_colors <- node_palette[round(contrast_norm * 99) + 1]

# Edge colors: light to dark grey based on validation strength
edge_weights <- igraph::E(network_graph)$weight
edge_palette <- grDevices::colorRampPalette(c("#d9d9d9", "#252525"))(100)
norm_weights <- (edge_weights - min(edge_weights)) / (max(edge_weights) - min(edge_weights))
edge_colors <- edge_palette[round(norm_weights * 99) + 1]

# Edge widths proportional to validation count
edge_widths <- edge_weights * 2

# Label colors: black for light nodes, white for dark nodes
label_colors <- ifelse(contrast_norm < 0.2, "black", "white")

# Set seed for reproducible layout (Fruchterman-Reingold is stochastic)
set.seed(42)
layout_fr <- igraph::layout_with_fr(network_graph)

pdf(file = output_figure_file, width = 12, height = 12)
plot(
  network_graph,
  edge.width = edge_widths,
  edge.color = edge_colors,
  vertex.size = 15,
  vertex.label.color = label_colors,
  vertex.label.cex = 0.8,
  vertex.color = vertex_colors,
  layout = layout_fr
)

# Add legend for node degree (connection count)
degree_range <- range(node_degree)
legend(
  "bottomleft",
  legend = c(
    paste0(degree_range[1], " connections"),
    paste0(round(mean(degree_range)), " connections"),
    paste0(degree_range[2], " connections")
  ),
  fill = node_palette[c(1, 50, 100)],
  title = "Node Degree",
  cex = 1.5,
  bty = "n"
)

# Add legend for edge weight (validation strength)
unique_weights <- sort(unique(edge_weights))
legend(
  "bottomright",
  legend = paste(unique_weights, "cohorts"),
  fill = edge_palette[round((unique_weights - min(edge_weights)) / (max(edge_weights) - min(edge_weights)) * 99) + 1],
  title = "Validation Strength",
  cex = 1.5,
  bty = "n"
)
dev.off()
