## SPDX-FileCopyrightText: 2026 Alkahest, Inc. <chronospd_data_access@alkahest.com>
##
## SPDX-License-Identifier: MIT

# ------------------------------------------------------------------------------
# Function: cluster_terms_simplify
# Description: Wrapper for simplifyEnrichment package to cluster terms based on gene association similarity
# Parameters:
#   term_ids - character vector of term IDs to cluster
#   TG - gene x term matrix (sparse or dense)
#   gene_ids - optional subset of genes (default: NULL = all)
#   method - clustering method (default: "binary_cut")
#   control - list passed to cluster_terms() (default: list())
#   min_set_size - minimum gene set size (default: 1)
#   sim_metric - similarity metric: "jaccard", "dice", "overlap", "kappa" (default: "jaccard")
#   kappa_transform - kappa transform: "clip", "rescale", "none" (default: "clip")
#   verbose - print progress (default: TRUE)
# Returns: list with clusters, similarity matrix, dropped terms, method info
# ------------------------------------------------------------------------------
cluster_terms_simplify <- function(
  term_ids,
  TG,
  gene_ids = NULL,
  method = "binary_cut",
  control = list(),
  min_set_size = 1,
  sim_metric = c("jaccard", "dice", "overlap", "kappa"),
  kappa_transform = c("clip", "rescale", "none"),
  verbose = TRUE
) {
  
  # ============================================================================
  # Internal helper function
  # ============================================================================
  
  # Compute term-term similarity matrix
  compute_similarity = function(gene_term_matrix, set_sizes, metric, kappa_transform) {
    intersection = as.matrix(Matrix::crossprod(gene_term_matrix))
    sizes = as.numeric(set_sizes)
    n_terms = length(sizes)
    
    # Replicate sizes as matrices
    size_a = matrix(sizes, nrow = n_terms, ncol = n_terms)
    size_b = t(size_a)
    
    similarity = switch(
      metric,
      "jaccard" = {
        denominator = size_a + size_b - intersection
        ifelse(denominator > 0, intersection / denominator, 0)
      },
      "dice" = {
        denominator = size_a + size_b
        ifelse(denominator > 0, (2 * intersection) / denominator, 0)
      },
      "overlap" = {
        denominator = pmin(size_a, size_b)
        ifelse(denominator > 0, intersection / denominator, 0)
      },
      "kappa" = {
        n_genes = nrow(gene_term_matrix)
        d = n_genes - (size_a + size_b - intersection)
        p_observed = (intersection + d) / n_genes
        p_expected = (size_a * size_b + (n_genes - size_a) * (n_genes - size_b)) / (n_genes * n_genes)
        result = (p_observed - p_expected) / (1 - p_expected)
        result[!is.finite(result)] = 0
        
        if (kappa_transform == "clip") {
          result[result < 0] = 0
        } else if (kappa_transform == "rescale") {
          result = (result + 1) / 2
        }
        result
      }
    )
    
    diag(similarity) = 1
    rownames(similarity) = colnames(gene_term_matrix)
    colnames(similarity) = colnames(gene_term_matrix)
    similarity
  }
  
  # ============================================================================
  # Validate inputs
  # ============================================================================
  
  if (!inherits(TG, c("matrix", "dgCMatrix", "dgTMatrix", "dgeMatrix"))) {
    stop("term clustering: TG must be a matrix")
  }
  
  if (is.null(rownames(TG)) || is.null(colnames(TG))) {
    stop("term clustering: TG must have rownames (entities) and colnames (terms)")
  }
  
  sim_metric = match.arg(sim_metric)
  kappa_transform = match.arg(kappa_transform)
  
  # Check if method is valid
  available_methods = try(simplifyEnrichment::all_clustering_methods(), silent = TRUE)
  if (!inherits(available_methods, "try-error") && !(method %in% available_methods)) {
    stop(sprintf("term clustering: cannot recognize method '%s'. Available: %s", method, paste(available_methods, collapse = ", ")))
  }
  
  term_ids = unique(as.character(term_ids))
  
  # ============================================================================
  # Orient matrix to genes x terms
  # ============================================================================
  
  in_cols = sum(term_ids %in% colnames(TG))
  in_rows = sum(term_ids %in% rownames(TG))
  if (in_cols == 0 && in_rows == 0) {
    stop("term clustering: none of the 'term_ids' present in TG rownames / colnames")
  }
  
  # Transpose if needed
  gene_term_matrix = if (in_cols >= in_rows) TG else t(TG)
  
  # ============================================================================
  # Trim to requested terms and genes
  # ============================================================================
  
  present_terms = intersect(term_ids, colnames(gene_term_matrix))
  dropped_terms = setdiff(term_ids, present_terms)
  
  if (!length(present_terms)) {
    warning("term clustering: no query terms present after matching; returning empty result")
    return(list(
      clusters = data.frame(term.id = character(), cluster = integer(), set_size = integer()),
      sim = matrix(numeric(0), 0, 0),
      dropped_terms = unique(dropped_terms),
      method = method,
      control = control,
      used_fallback = NA,
      sim_metric = sim_metric
    ))
  }
  
  if (!is.null(gene_ids)) {
    present_genes = intersect(as.character(gene_ids), rownames(gene_term_matrix))
    if (!length(present_genes)) {
      warning("term clustering: no query genes present after matching; empty result.")
      return(list(
        clusters = data.frame(term.id = character(), cluster = integer(), set_size = integer()),
        sim = matrix(numeric(0), 0, 0),
        dropped_terms = unique(c(dropped_terms, present_terms)),
        method = method,
        control = control,
        used_fallback = NA,
        sim_metric = sim_metric
      ))
    }
    gene_term_matrix = gene_term_matrix[present_genes, present_terms, drop = FALSE]
  } else {
    gene_term_matrix = gene_term_matrix[, present_terms, drop = FALSE]
  }
  
  # ============================================================================
  # Binarize matrix
  # ============================================================================
  
  if (inherits(gene_term_matrix, "dgCMatrix")) {
    gene_term_matrix = methods::as(gene_term_matrix, "dgCMatrix")
    if (length(gene_term_matrix@x)) gene_term_matrix@x[] = 1
    gene_term_matrix = Matrix::drop0(gene_term_matrix)
  } else if (is.matrix(gene_term_matrix)) {
    gene_term_matrix = Matrix::Matrix((gene_term_matrix != 0) + 0, sparse = TRUE)
  } else {
    gene_term_matrix = Matrix::Matrix(as.matrix(gene_term_matrix) != 0, sparse = TRUE)
  }
  
  # ============================================================================
  # Filter by set size
  # ============================================================================
  
  set_sizes = Matrix::colSums(gene_term_matrix)
  names(set_sizes) = colnames(gene_term_matrix)
  keep_terms = set_sizes >= as.integer(min_set_size)
  
  if (!any(keep_terms)) {
    warning("term clustering: no terms remain after min_set_size filtering; returning empty result.")
    return(list(
      clusters = data.frame(term.id = character(), cluster = integer(), set_size = integer()),
      sim = matrix(numeric(0), 0, 0),
      dropped_terms = unique(c(dropped_terms, colnames(gene_term_matrix))),
      method = method,
      control = control,
      used_fallback = NA,
      sim_metric = sim_metric
    ))
  }
  
  if (!all(keep_terms)) {
    dropped_terms = unique(c(dropped_terms, colnames(gene_term_matrix)[!keep_terms]))
    gene_term_matrix = gene_term_matrix[, keep_terms, drop = FALSE]
    set_sizes = set_sizes[colnames(gene_term_matrix)]
  }
  
  term_names = colnames(gene_term_matrix)
  if (!length(term_names)) {
    warning("term clustering: no terms left after filtering; returning empty result.")
    return(list(
      clusters = data.frame(term.id = character(), cluster = integer(), set_size = integer()),
      sim = matrix(numeric(0), 0, 0),
      dropped_terms = unique(dropped_terms),
      method = method,
      control = control,
      used_fallback = NA,
      sim_metric = sim_metric
    ))
  }
  
  # ============================================================================
  # Handle single term case
  # ============================================================================
  
  if (ncol(gene_term_matrix) == 1L) {
    similarity = matrix(1, 1, 1, dimnames = list(term_names, term_names))
    output = data.frame(
      term.id = term_names,
      cluster = 1L,
      set_size = as.integer(set_sizes[term_names]),
      row.names = NULL,
      check.names = FALSE
    )
    return(list(
      clusters = output,
      sim = similarity,
      dropped_terms = unique(dropped_terms),
      method = method,
      control = control,
      used_fallback = FALSE,
      sim_metric = sim_metric
    ))
  }
  
  # ============================================================================
  # Compute similarity matrix
  # ============================================================================
  
  similarity = compute_similarity(gene_term_matrix, set_sizes, sim_metric, kappa_transform)
  
  # ============================================================================
  # Try clustering with simplifyEnrichment
  # ============================================================================
  
  used_fallback = FALSE
  clusters = NULL
  error_message = NULL
  
  clusters = tryCatch(
    simplifyEnrichment::cluster_terms(similarity, method = method, control = control, verbose = verbose),
    error = function(e) {
      error_message <<- conditionMessage(e)
      NULL
    }
  )
  
  # ============================================================================
  # Fallback to connected components if clustering failed
  # ============================================================================
  
  if (is.null(clusters) || length(clusters) != ncol(similarity)) {
    used_fallback = TRUE
    control_str = tryCatch(utils::capture.output(str(control)), error = function(e) "<unprintable>")
    warning(sprintf(
      "cluster_terms() failed (method='%s'). Falling back to simple connected components on sim>0.\nError: %s\ncontrol=:\n%s",
      method,
      if (is.null(error_message)) "unknown error" else error_message,
      paste(control_str, collapse = "\n")
    ))
    
    if (!requireNamespace("igraph", quietly = TRUE)) {
      stop("term clustering: need igraph package!")
    }
    
    graph = igraph::graph_from_adjacency_matrix(similarity > 0, mode = "undirected", diag = FALSE)
    membership = igraph::components(graph)$membership
    clusters = membership[match(term_names, names(membership))]
    names(clusters) = term_names
  } else {
    # Align names
    if (!is.null(names(clusters)) && setequal(names(clusters), term_names)) {
      clusters = clusters[term_names]
    } else {
      names(clusters) = term_names
    }
  }
  
  # ============================================================================
  # Normalize cluster labels to 1..K
  # ============================================================================
  
  if (is.factor(clusters)) clusters = as.character(clusters)
  clusters = if (!is.numeric(clusters)) {
    as.integer(factor(clusters, levels = unique(clusters)))
  } else {
    as.integer(factor(clusters, levels = unique(clusters)))
  }
  names(clusters) = term_names
  
  cluster_df = data.frame(
    term.id = term_names,
    cluster = clusters,
    set_size = as.integer(set_sizes[term_names]),
    stringsAsFactors = FALSE,
    row.names = NULL,
    check.names = FALSE
  )
  cluster_df = cluster_df[order(cluster_df$cluster, -cluster_df$set_size, cluster_df$term.id), , drop = FALSE]
  
  # ============================================================================
  # Output
  # ============================================================================
  
  list(
    clusters = cluster_df,
    sim = similarity,
    dropped_terms = unique(dropped_terms),
    method = method,
    control = control,
    used_fallback = used_fallback,
    sim_metric = sim_metric,
    kappa_transform = kappa_transform
  )
}
