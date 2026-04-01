## SPDX-FileCopyrightText: 2026 Alkahest, Inc. <chronospd_data_access@alkahest.com>
##
## SPDX-License-Identifier: MIT

# ------------------------------------------------------------------------------
# Function: summarize_PA_matrix
# Description: Extract PA run results to matrix formats for easier manipulation
# Parameters:
#   PA - PA run result object
# Returns: list with terms data frame and matrices (n.hits, p.adj, p.raw, fc)
# ------------------------------------------------------------------------------
summarize_PA_matrix <- function(PA) {
  
  # Get all terms from deepest set to build framework for matrix
  n_depths = length(PA$tested_entities)
  last_depth_df = do.call(rbind, list(PA$GO[[n_depths]], PA$Reactome[[n_depths]], PA$KEGG[[n_depths]]))
  all_terms = last_depth_df[, c("term.id", "AnnotDb", "Description", "n.background")]
  term_ids = all_terms$term.id
  
  # Prepare frames
  column_names = sapply(PA$tested_entities, nrow)
  term_list = vector("list", n_depths)
  
  for (k in seq_len(n_depths)) {
    depth_df = do.call(rbind, list(PA$GO[[k]], PA$Reactome[[k]], PA$KEGG[[k]]))
    values = depth_df[, c("term.id", "n.hit", "fdr", "p.value", "FoldChange")]
    column_label = as.character(nrow(PA$tested_entities[[k]]))
    
    term_list[[k]] = list(
      n.hit = setNames(values$n.hit, values$term.id),
      fdr = setNames(values$fdr, values$term.id),
      p.value = setNames(values$p.value, values$term.id),
      FoldChange = setNames(values$FoldChange, values$term.id),
      colname = column_label
    )
  }
  
  # Initialize matrices with NA placeholders
  matrix_init = matrix(
    NA_real_,
    nrow = length(term_ids),
    ncol = n_depths,
    dimnames = list(term_ids, sapply(term_list, `[[`, "colname"))
  )
  
  result = list(
    terms = all_terms,
    n.hits = matrix_init,
    p.adj = matrix_init,
    p.raw = matrix_init,
    fc = matrix_init
  )
  
  # Fill matrices with values
  for (k in seq_len(n_depths)) {
    result$n.hits[, k] = term_list[[k]]$n.hit[term_ids]
    result$p.adj[, k] = term_list[[k]]$fdr[term_ids]
    result$p.raw[, k] = term_list[[k]]$p.value[term_ids]
    result$fc[, k] = term_list[[k]]$FoldChange[term_ids]
  }
  
  return(result)
}
