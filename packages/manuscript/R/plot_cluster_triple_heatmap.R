## SPDX-FileCopyrightText: 2026 Alkahest, Inc. <chronospd_data_access@alkahest.com>
##
## SPDX-License-Identifier: MIT

# ------------------------------------------------------------------------------
# Function: plot_cluster_triple_heatmap
# Description: Generate triple heatmap showing pathway clusters across conditions (All, preDx, postDx)
# Parameters:
#   sel_terms - term IDs used for clustering
#   clustering_result - list with clustering_result$clusters (term.id, cluster)
#   pa_matrices_synced - list with named PA matrix objects
#   drange - SEPA depth range as column indices
#   fdr_cut - FDR cutoff for visualization (default: 0.05)
#   cluster_names - optional pretty names per cluster id (default: NULL)
#   append_n - append "(n = XXX terms)" to cluster labels (default: TRUE)
#   show_row_names - show cluster names (default: TRUE)
#   rect_gp_col - rectangle border color (default: "grey20")
#   row_title_rot - row title rotation (default: 0)
#   titles - panel titles (default: c("All Samples","preDx","postDx"))
# Returns: invisible(NULL) - draws heatmap
# ------------------------------------------------------------------------------
plot_cluster_triple_heatmap <- function(
  sel_terms,
  clustering_result,
  pa_matrices_synced,
  drange,
  fdr_cut = 0.05,
  cluster_names = NULL,
  append_n = TRUE,
  show_row_names = TRUE,
  rect_gp_col = "grey20",
  row_title_rot = 0,
  titles = c("All Samples", "preDx", "postDx")
) {
  
  # Check required packages
  stopifnot(
    requireNamespace("ComplexHeatmap", quietly = TRUE),
    requireNamespace("circlize", quietly = TRUE),
    requireNamespace("grid", quietly = TRUE)
  )
  
  # ============================================================================
  # Internal helper functions
  # ============================================================================
  
  # Pull FDR matrix from PA condition within depth range, align to term_ids
  safe_get_fdr = function(pa_slot, term_ids, drange) {
    fdr_matrix = pa_slot$p.adj
    term_ids_slot = pa_slot$terms$term.id
    match_indices = match(term_ids, term_ids_slot)
    
    output_matrix = matrix(
      NA_real_,
      nrow = length(term_ids),
      ncol = length(drange)
    )
    rownames(output_matrix) = term_ids
    colnames(output_matrix) = colnames(fdr_matrix)[drange]
    
    ok_indices = which(!is.na(match_indices))
    if (length(ok_indices)) {
      output_matrix[ok_indices, ] = as.matrix(fdr_matrix[match_indices[ok_indices], drange, drop = FALSE])
    }
    output_matrix
  }
  
  # Minimum ignoring non-finite
  safe_min = function(vec) {
    vec = vec[is.finite(vec)]
    if (!length(vec)) NA_real_ else min(vec)
  }
  
  # Aggregate up/down FDR to cluster level for one panel
  aggregate_panel = function(fdr_up, fdr_down, cluster_ids, terms_by_cluster, cluster_labels, fdr_cut) {
    
    n_clusters = length(cluster_ids)
    threshold = -log10(fdr_cut)
    
    # Initialize output matrix
    output_matrix = matrix(
      0,
      nrow = n_clusters,
      ncol = ncol(fdr_up),
      dimnames = list(as.vector(cluster_labels[as.character(cluster_ids)]), colnames(fdr_up))
    )
    
    # For each cluster
    for (k in seq_along(cluster_ids)) {
      cluster_id = cluster_ids[k]
      cluster_terms = terms_by_cluster[[as.character(cluster_id)]]
      if (!length(cluster_terms)) next
      
      # Extract up/down FDR for terms in this cluster
      up_fdr = fdr_up[cluster_terms, , drop = FALSE]
      down_fdr = fdr_down[cluster_terms, , drop = FALSE]
      
      # Best (minimum) FDR within cluster for each depth
      up_min = apply(up_fdr, 2, safe_min)
      down_min = apply(down_fdr, 2, safe_min)
      
      # Convert to signed values (-log10)
      up_values = -log10(up_min)
      down_values = log10(down_min)
      
      # Handle NA
      up_values[!is.finite(up_values)] = 0
      down_values[!is.finite(down_values)] = 0
      
      # Merge: pick stronger absolute direction
      merged_values = ifelse(abs(up_values) >= abs(down_values), up_values, down_values)
      
      # Set non-significant to 0
      merged_values[is.na(merged_values) | abs(merged_values) < threshold] = 0
      output_matrix[k, ] = merged_values
    }
    output_matrix
  }
  
  # Convert to numeric matrix for ComplexHeatmap
  to_numeric_matrix = function(x) {
    x = as.matrix(x)
    storage.mode(x) = "numeric"
    x[is.na(x)] = 0
    x
  }
  
  # Wrapper for ComplexHeatmap::Heatmap with settings
  make_heatmap = function(matrix, title, show_dendrogram, legend_name, with_legend = FALSE) {
    ComplexHeatmap::Heatmap(
      matrix,
      name = legend_name,
      col = color_function,
      column_title = title,
      cluster_rows = isTRUE(show_dendrogram),
      cluster_row_slices = FALSE,
      row_split = row_split_factor,
      row_title = NULL,
      show_row_dend = isTRUE(show_dendrogram),
      cluster_columns = FALSE,
      show_row_names = show_row_names,
      show_column_names = TRUE,
      row_names_gp = grid::gpar(fontsize = 8),
      column_names_gp = grid::gpar(fontsize = 8),
      rect_gp = grid::gpar(col = rect_gp_col, lwd = 1),
      row_gap = grid::unit(3, "mm"),
      width = grid::unit(ncol(matrix) * 0.25, "cm"),
      show_heatmap_legend = with_legend,
      heatmap_legend_param = if (with_legend) {
        list(title = "-log10 FDR", direction = "horizontal")
      } else {
        list()
      }
    )
  }
  
  # ============================================================================
  # Prepare clusters and terms
  # ============================================================================
  
  sel_terms = unique(as.character(sel_terms))
  cluster_df = clustering_result$clusters[, c("term.id", "cluster")]
  cluster_df$term.id = as.character(cluster_df$term.id)
  cluster_df$cluster = as.integer(cluster_df$cluster)
  
  term_ids = intersect(sel_terms, cluster_df$term.id)
  if (!length(term_ids)) {
    stop("No overlap between 'sel_terms' and clustering_result$clusters$term.id.")
  }
  
  # ============================================================================
  # Collect term-level FDR
  # ============================================================================
  
  fdr_all_up = safe_get_fdr(pa_matrices_synced[[PA_CONDITION_ALL_UP]], term_ids, drange)
  fdr_all_down = safe_get_fdr(pa_matrices_synced[[PA_CONDITION_ALL_DOWN]], term_ids, drange)
  fdr_pre_up = safe_get_fdr(pa_matrices_synced[[PA_CONDITION_PRE_UP]], term_ids, drange)
  fdr_pre_down = safe_get_fdr(pa_matrices_synced[[PA_CONDITION_PRE_DOWN]], term_ids, drange)
  fdr_post_up = safe_get_fdr(pa_matrices_synced[[PA_CONDITION_POST_UP]], term_ids, drange)
  fdr_post_down = safe_get_fdr(pa_matrices_synced[[PA_CONDITION_POST_DOWN]], term_ids, drange)
  
  if (!any(is.finite(fdr_all_up)) &&
      !any(is.finite(fdr_all_down)) &&
      !any(is.finite(fdr_pre_up)) &&
      !any(is.finite(fdr_pre_down)) &&
      !any(is.finite(fdr_post_up)) &&
      !any(is.finite(fdr_post_down))) {
    stop("heatmap: no finite FDR values in any PA panel for the selected terms.")
  }
  
  # ============================================================================
  # Sync clusters to selected terms
  # ============================================================================
  
  cluster_df = cluster_df[cluster_df$term.id %in% term_ids, , drop = FALSE]
  if (!nrow(cluster_df)) {
    stop("After matching to selected terms, no clustered terms remain.")
  }
  
  cluster_ids = sort(unique(cluster_df$cluster))
  terms_by_cluster = split(cluster_df$term.id, cluster_df$cluster)
  
  # ============================================================================
  # Build cluster labels
  # ============================================================================
  
  cluster_labels = setNames(paste0("Cluster ", cluster_ids), cluster_ids)
  
  if (!is.null(cluster_names)) {
    if (is.data.frame(cluster_names) && all(c("cluster", "label") %in% names(cluster_names))) {
      temp_labels = setNames(as.character(cluster_names$label), as.character(cluster_names$cluster))
      cluster_labels[names(temp_labels)] = temp_labels
    } else if (!is.null(names(cluster_names))) {
      cluster_labels[names(cluster_names)] = as.character(cluster_names)
    }
  }
  
  if (append_n) {
    for (cluster_id in cluster_ids) {
      n_terms_in_cluster = length(terms_by_cluster[[as.character(cluster_id)]])
      cluster_labels[as.character(cluster_id)] = paste0(cluster_labels[as.character(cluster_id)], " (n=", n_terms_in_cluster, ")")
    }
  }
  
  # ============================================================================
  # Aggregate to cluster level
  # ============================================================================
  
  matrix_all = to_numeric_matrix(aggregate_panel(fdr_all_up, fdr_all_down, cluster_ids, terms_by_cluster, cluster_labels, fdr_cut))
  matrix_pre = to_numeric_matrix(aggregate_panel(fdr_pre_up, fdr_pre_down, cluster_ids, terms_by_cluster, cluster_labels, fdr_cut))
  matrix_post = to_numeric_matrix(aggregate_panel(fdr_post_up, fdr_post_down, cluster_ids, terms_by_cluster, cluster_labels, fdr_cut))
  
  # ============================================================================
  # Classify by directionality
  # ============================================================================
  
  any_positive = (rowSums(matrix_all > 0) + rowSums(matrix_pre > 0) + rowSums(matrix_post > 0)) > 0
  any_negative = (rowSums(matrix_all < 0) + rowSums(matrix_pre < 0) + rowSums(matrix_post < 0)) > 0
  
  directionality = ifelse(
    any_positive & any_negative,
    "BI-DIRECTIONAL",
    ifelse(any_positive, "UPREGULATED", ifelse(any_negative, "DOWNREGULATED", "NEUTRAL"))
  )
  
  desired_order = c("UPREGULATED", "DOWNREGULATED", "BI-DIRECTIONAL", "NEUTRAL")
  present_levels = desired_order[desired_order %in% unique(directionality)]
  row_split_factor = factor(directionality, levels = present_levels)
  
  # ============================================================================
  # Build slice labels with term counts
  # ============================================================================
  
  n_terms_per_slice = vapply(present_levels, function(slice_level) {
    cluster_ids_in_slice = cluster_ids[row_split_factor == slice_level]
    if (!length(cluster_ids_in_slice)) return(0L)
    all_terms_in_slice = unlist(terms_by_cluster[as.character(cluster_ids_in_slice)], use.names = FALSE)
    length(unique(all_terms_in_slice))
  }, integer(1))
  
  slice_labels = paste0(present_levels, "\n(n=", n_terms_per_slice, ")")
  
  # Left annotation with slice labels
  left_annotation = ComplexHeatmap::rowAnnotation(
    groups = ComplexHeatmap::anno_block(
      labels = slice_labels,
      labels_rot = row_title_rot,
      labels_gp = grid::gpar(fontsize = 10),
      gp = grid::gpar(fill = NA, col = NA)
    )
  )
  
  # ============================================================================
  # Color scale
  # ============================================================================
  
  color_function = circlize::colorRamp2(
    c(-6, log10(0.05), log10(0.0501), -log10(0.0501), -log10(0.05), 6),
    c("navyblue", "turquoise", "gray90", "gray90", "yellow", "darkred")
  )
  
  # ============================================================================
  # Build heatmap panels
  # ============================================================================
  
  heatmap_all = make_heatmap(matrix_all, titles[1], show_dendrogram = TRUE, legend_name = "FDR_all", with_legend = TRUE)
  heatmap_pre = make_heatmap(matrix_pre, titles[2], show_dendrogram = FALSE, legend_name = "FDR_pre", with_legend = FALSE)
  heatmap_post = make_heatmap(matrix_post, titles[3], show_dendrogram = FALSE, legend_name = "FDR_post", with_legend = FALSE)
  
  # ============================================================================
  # Draw composite heatmap
  # ============================================================================
  
  ComplexHeatmap::draw(
    left_annotation + heatmap_all + heatmap_pre + heatmap_post,
    merge_legend = TRUE,
    heatmap_legend_side = "top",
    ht_gap = grid::unit(5, "mm"),
    padding = grid::unit(c(2, 2, 2, 20), "mm")
  )
  
  invisible(NULL)
}
