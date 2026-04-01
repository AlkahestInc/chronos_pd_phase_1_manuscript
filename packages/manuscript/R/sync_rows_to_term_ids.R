## SPDX-FileCopyrightText: 2026 Alkahest, Inc. <chronospd_data_access@alkahest.com>
##
## SPDX-License-Identifier: MIT

# ------------------------------------------------------------------------------
# Function: sync_rows_to_term_ids
# Description: Sync PA matrices to a global set of terms (same row = same term; missing is NA)
# Parameters:
#   pa_matrix - PA matrix result from summarize_PA_matrix()
#   term_ids - character vector of term IDs to sync to
#   matrix_names - names of matrices to sync (default: c("n.hits", "p.raw", "p.adj", "fc"))
# Returns: synced PA matrix with rows aligned to term_ids
# ------------------------------------------------------------------------------
sync_rows_to_term_ids <- function(pa_matrix, term_ids, matrix_names = c("n.hits", "p.raw", "p.adj", "fc")) {
  stopifnot("terms" %in% names(pa_matrix), "term.id" %in% colnames(pa_matrix$terms))
  
  match_indices = match(term_ids, pa_matrix$terms$term.id)
  matched_rows = !is.na(match_indices)
  
  # Expand terms
  template = pa_matrix$terms[0, , drop = FALSE]
  full_terms = as.data.frame(matrix(NA, nrow = length(term_ids), ncol = ncol(template)))
  colnames(full_terms) = colnames(template)
  full_terms$term.id = term_ids
  full_terms[matched_rows, ] = pa_matrix$terms[match_indices[matched_rows], , drop = FALSE]
  rownames(full_terms) = term_ids
  
  # Initialize output
  output = list(terms = full_terms)
  
  # Expand matrices
  for (matrix_name in matrix_names) {
    if (!is.null(pa_matrix[[matrix_name]])) {
      expanded_matrix = matrix(NA, nrow = length(term_ids), ncol = ncol(pa_matrix[[matrix_name]]))
      colnames(expanded_matrix) = colnames(pa_matrix[[matrix_name]])
      rownames(expanded_matrix) = term_ids
      expanded_matrix[matched_rows, ] = pa_matrix[[matrix_name]][match_indices[matched_rows], , drop = FALSE]
      output[[matrix_name]] = expanded_matrix
    }
  }
  
  return(output)
}
