## SPDX-FileCopyrightText: 2026 Alkahest, Inc. <chronospd_data_access@alkahest.com>
##
## SPDX-License-Identifier: MIT

# ------------------------------------------------------------------------------
# Function: summarize_significant_terms
# Description: Extract compressed SEPA summary - one row per pathway
# Parameters:
#   PA - PA run result object
#   pcut - p-value cutoff (default: 0.05)
#   ncut - minimum number of hits (default: 3)
#   fccut - fold change cutoff (default: 0)
#   scut - required consecutive SEPA steps (default: 3)
#   sizecut - term size range c(min, max) (default: c(0, 500))
#   useFDR - use FDR instead of raw p-value (default: TRUE)
# Returns: data frame with significant terms summary
# ------------------------------------------------------------------------------
summarize_significant_terms <- function(PA, pcut = 0.05, ncut = 3, fccut = 0,
                                       scut = 3, sizecut = c(0, 500), useFDR = TRUE) {
  
  # Convert to matrix format
  pa_matrices = summarize_PA_matrix(PA)
  n_depths = length(PA$tested_entities)
  
  # Choose p-value or FDR filter
  p_values = if (useFDR) pa_matrices$p.adj else pa_matrices$p.raw
  
  # Binary mask for passing criteria
  passed = p_values <= pcut & pa_matrices$n.hits >= ncut & pa_matrices$fc >= fccut
  passed[is.na(passed)] = FALSE
  
  # SEPA core filter: consecutive TRUE across depths in a window
  keep_terms = apply(passed, 1, function(x) {
    rle_result = rle(x)
    any(rle_result$lengths[rle_result$values] >= scut)
  })
  
  size_ok = pa_matrices$terms$n.background >= sizecut[1] & pa_matrices$terms$n.background <= sizecut[2]
  passed_term_indices = which(keep_terms & size_ok)
  
  # Return empty if no terms passed
  if (length(passed_term_indices) == 0) {
    message("No terms passed filtering thresholds.")
    return(data.frame())
  }
  
  # Add annotation list
  annotation_list = Map(function(go, reactome, kegg) {
    combined_df = do.call(rbind, list(go, reactome, kegg))
    combined_df$term.id = as.character(combined_df$term.id)
    combined_df
  }, PA$GO, PA$Reactome, PA$KEGG)
  
  # Column labels (SEPA depth)
  column_ids = colnames(p_values)
  selected_term_ids = rownames(p_values)[passed_term_indices]
  
  # Assemble summary
  result_list = lapply(selected_term_ids, function(term_id) {
    pvals_raw = pa_matrices$p.raw[term_id, ]
    pvals_adj = pa_matrices$p.adj[term_id, ]
    n_hits = pa_matrices$n.hits[term_id, ]
    fold_changes = pa_matrices$fc[term_id, ]
    
    valid_depths = which(!is.na(pvals_raw))
    if (length(valid_depths) == 0) return(NULL)
    
    first_pass_depth = which(passed[term_id, ])[1]
    
    p_used = if (useFDR) pvals_adj else pvals_raw
    p_used[!is.finite(p_used)] = Inf
    best_pval_depth = which.min(p_used)
    
    max_depth_index = max(valid_depths)
    
    annotation_row = annotation_list[[max_depth_index]]
    annotation_row = annotation_row[annotation_row$term.id == term_id, , drop = FALSE]
    
    entrez_ids = if ("ENTREZID" %in% names(annotation_row)) {
      unique(unlist(strsplit(annotation_row$ENTREZID, "/", fixed = TRUE)))
    } else NA
    
    symbol_ids = if ("SYMBOL" %in% names(annotation_row)) {
      unique(unlist(strsplit(annotation_row$SYMBOL, "/", fixed = TRUE)))
    } else NA
    
    # Match metadata using term.id as key
    term_metadata = pa_matrices$terms[pa_matrices$terms$term.id == term_id, , drop = FALSE]
    
    data.frame(
      term.id = term_id,
      AnnotDb = term_metadata$AnnotDb,
      Description = term_metadata$Description,
      n.background = term_metadata$n.background,
      first_depth_pass = column_ids[first_pass_depth],
      best_pval_depth = column_ids[best_pval_depth],
      p.raw_at_best = pvals_raw[best_pval_depth],
      p.adj_at_best = pvals_adj[best_pval_depth],
      n.hit_at_best = n_hits[best_pval_depth],
      fc_at_best = fold_changes[best_pval_depth],
      depth_span = length(valid_depths),
      ENTREZID = paste(sort(na.omit(entrez_ids)), collapse = "/"),
      SYMBOL = paste(sort(na.omit(symbol_ids)), collapse = "/"),
      stringsAsFactors = FALSE
    )
  })
  
  summary_df = do.call(rbind, result_list)
  rownames(summary_df) = NULL
  return(summary_df)
}
