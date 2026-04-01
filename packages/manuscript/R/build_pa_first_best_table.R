## SPDX-FileCopyrightText: 2026 Alkahest, Inc. <chronospd_data_access@alkahest.com>
##
## SPDX-License-Identifier: MIT

# ------------------------------------------------------------------------------
# Function: build_pa_first_best_table
# Description: Extended SEPA compressor for supplementary tables - can report 2 rows per pathway
#              (first appearance, most significant)
# Parameters:
#   PA - PA run result object
#   pcut - p-value cutoff (default: 0.05)
#   ncut - minimum hits (default: 2)
#   fccut - fold change cutoff (default: 0)
#   scut - consecutive steps (default: 3)
#   sizecut - term size range c(min, max) (default: c(0, Inf))
#   drange - SEPA depth range (default: NULL = all)
#   useFDR - use FDR (default: TRUE)
#   annotdbs - annotation databases to include (default: NULL = all)
#   prefer_symbol - prefer SYMBOL over ENTREZID for hits (default: TRUE)
# Returns: data frame with first/best rows per pathway
# ------------------------------------------------------------------------------
build_pa_first_best_table <- function(PA, pcut = 0.05, ncut = 2, fccut = 0,
                                     scut = 3, sizecut = c(0, Inf), drange = NULL, useFDR = TRUE,
                                     annotdbs = NULL, prefer_symbol = TRUE) {
  
  # ============================================================================
  # Internal helper functions
  # ============================================================================
  
  # Find regions of consecutive TRUE values
  find_regions = function(x, min_len) {
    rle_result = rle(as.logical(x))
    ends = cumsum(rle_result$lengths)
    starts = ends - rle_result$lengths + 1
    ok_indices = which(rle_result$values & rle_result$lengths >= min_len)
    if (!length(ok_indices)) return(list(starts = integer(0), ends = integer(0)))
    list(starts = starts[ok_indices], ends = ends[ok_indices])
  }
  
  # Translate drange-relative positions to topn labels
  rel_to_topn = function(rel_idx_vec) column_labels[depth_columns[rel_idx_vec]]
  
  # Build one row for output table
  make_one_row = function(rel_idx, abs_idx, is_first_flag, is_best_flag) {
    region_index = which(rel_idx >= regions$starts & rel_idx <= regions$ends)
    if (length(region_index)) {
      region_start = regions$starts[region_index[1]]
      region_end = regions$ends[region_index[1]]
      region_length = region_end - region_start + 1L
      region_start_topn = rel_to_topn(region_start)
      region_end_topn = rel_to_topn(region_end)
    } else {
      region_start_topn = NA_character_
      region_end_topn = NA_character_
      region_length = NA_integer_
    }
    
    p_raw = pa_matrices$p.raw[current_term, abs_idx]
    p_adj = pa_matrices$p.adj[current_term, abs_idx]
    n_hit = pa_matrices$n.hits[current_term, abs_idx]
    fc = pa_matrices$fc[current_term, abs_idx]
    topn = column_labels[abs_idx]
    
    # Hits at selected depth
    annotation_depth = annotations_by_depth[[abs_idx]]
    if (!nrow(annotation_depth)) {
      hits_entrez = character(0)
      hits_symbol = character(0)
    } else {
      annotation_subset = annotation_depth[annotation_depth$term.id == current_term, , drop = FALSE]
      if (!nrow(annotation_subset)) {
        hits_entrez = character(0)
        hits_symbol = character(0)
      } else {
        entrez_string = as.character(annotation_subset$ENTREZID[1])
        symbol_string = as.character(annotation_subset$SYMBOL[1])
        hits_entrez = if (is.na(entrez_string) || !nzchar(entrez_string)) {
          character(0)
        } else {
          strsplit(entrez_string, "/", fixed = TRUE)[[1]]
        }
        hits_symbol = if (is.na(symbol_string) || !nzchar(symbol_string)) {
          character(0)
        } else {
          strsplit(symbol_string, "/", fixed = TRUE)[[1]]
        }
        
        tested_entities = PA$tested_entities[[abs_idx]]
        if (!is.data.frame(tested_entities) || !"ENTREZID" %in% names(tested_entities)) {
          stop(sprintf("SEPA compressor error: no ENTREZID column in tested_entities[[%d]]", abs_idx))
        }
        tested_entrez = unique(as.character(na.omit(tested_entities$ENTREZID)))
        hits_entrez = intersect(hits_entrez, tested_entrez)
        
        if (prefer_symbol && !length(hits_symbol) && length(hits_entrez)) {
          if (requireNamespace("AnnotationDbi", quietly = TRUE) &&
              requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
            symbol_map = try(
              AnnotationDbi::select(
                org.Hs.eg.db::org.Hs.eg.db,
                keys = hits_entrez,
                keytype = "ENTREZID",
                columns = "SYMBOL"
              ),
              silent = TRUE
            )
            if (!inherits(symbol_map, "try-error") && nrow(symbol_map)) {
              hits_symbol = symbol_map$SYMBOL[match(hits_entrez, symbol_map$ENTREZID)]
              hits_symbol = unique(hits_symbol[!is.na(hits_symbol) & nzchar(hits_symbol)])
            }
          }
        }
      }
    }
    
    # Output row
    data.frame(
      term.id = current_term,
      AnnotDb = term_metadata$AnnotDb,
      Description = term_metadata$Description,
      n.background = term_metadata$n.background,
      is_first = is_first_flag,
      is_best = is_best_flag,
      topn = topn,
      p.raw = p_raw,
      p.adj = p_adj,
      n.hit = n_hit,
      fc = fc,
      region_start_topn = region_start_topn,
      region_end_topn = region_end_topn,
      region_len = region_length,
      hits.SYMBOL = if (length(hits_symbol)) paste(hits_symbol, collapse = "/") else NA_character_,
      hits.ENTREZID = if (length(hits_entrez)) paste(hits_entrez, collapse = "/") else NA_character_,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  }
  
  # ============================================================================
  # Sanity checks
  # ============================================================================
  
  stopifnot(is.list(PA))
  needed_elements = c("GO", "Reactome", "KEGG", "tested_entities", "dataset_background", "id_type")
  missing_elements = setdiff(needed_elements, names(PA))
  if (length(missing_elements)) {
    stop(sprintf("PA missing: %s", paste(missing_elements, collapse = ", ")))
  }
  
  if (!is.list(PA$GO) || !is.list(PA$Reactome) || !is.list(PA$KEGG)) {
    stop("SEPA compressor error: PA$GO / PA$Reactome / PA$KEGG must be lists with per-depth results")
  }
  
  if (!is.list(PA$tested_entities)) {
    stop("SEPA compressor error: PA$tested_entities must be a list")
  }
  
  n_depths = length(PA$tested_entities)
  if (!all(lengths(list(PA$GO, PA$Reactome, PA$KEGG)) == n_depths)) {
    stop("SEPA compressor error: lengths of PA$ lists must match")
  }
  
  # Switch to matrix format
  pa_matrices = summarize_PA_matrix(PA)
  
  needed_elements = c("terms", "n.hits", "p.raw", "p.adj", "fc")
  missing_elements = setdiff(needed_elements, names(pa_matrices))
  if (length(missing_elements)) {
    stop(sprintf("SEPA compressor error: summarize_PA_matrix(PA) missing: %s", paste(missing_elements, collapse = ", ")))
  }
  
  same_rows = all(
    rownames(pa_matrices$n.hits) == rownames(pa_matrices$p.raw),
    rownames(pa_matrices$p.raw) == rownames(pa_matrices$p.adj),
    rownames(pa_matrices$p.adj) == rownames(pa_matrices$fc)
  )
  same_cols = all(
    colnames(pa_matrices$n.hits) == colnames(pa_matrices$p.raw),
    colnames(pa_matrices$p.raw) == colnames(pa_matrices$p.adj),
    colnames(pa_matrices$p.adj) == colnames(pa_matrices$fc)
  )
  if (!same_rows || !same_cols) {
    stop("SEPA compressor error: mismatch in row/column names of PA matrices")
  }
  
  term_ids = rownames(pa_matrices$p.raw)
  column_labels = colnames(pa_matrices$p.raw)
  if (is.null(term_ids) || is.null(column_labels)) {
    stop("SEPA compressor error: PA rownames != term ids / colnames != depth labels")
  }
  
  # ============================================================================
  # Filter to SEPA depth range
  # ============================================================================
  
  if (is.null(drange)) {
    depth_columns = seq_len(ncol(pa_matrices$p.raw))
  } else if (is.numeric(drange)) {
    if (any(drange < 1 | drange > ncol(pa_matrices$p.raw))) {
      stop("SEPA compressor error: drange is not within SEPA depth range")
    }
    depth_columns = as.integer(drange)
  } else if (is.character(drange)) {
    match_indices = match(drange, column_labels)
    if (any(is.na(match_indices))) {
      stop(sprintf("SEPA compressor error: drange labels not found: %s", paste(drange[is.na(match_indices)], collapse = ", ")))
    }
    depth_columns = as.integer(match_indices)
  } else {
    stop("SEPA compressor error: drange must be NULL / numeric / character")
  }
  
  # ============================================================================
  # Filter to annotation databases
  # ============================================================================
  
  if (!is.null(annotdbs)) {
    if (!"AnnotDb" %in% names(pa_matrices$terms)) {
      stop("SEPA compressor error: 'AnnotDb' column not found in PA terms df")
    }
    keep_rows = pa_matrices$terms$AnnotDb %in% annotdbs
  } else {
    keep_rows = rep(TRUE, nrow(pa_matrices$terms))
  }
  
  # ============================================================================
  # Term size filter
  # ============================================================================
  
  needed_columns = c("n.background", "term.id", "AnnotDb", "Description")
  if (!all(needed_columns %in% names(pa_matrices$terms))) {
    stop("SEPA compressor error: terms df must have: term.id, AnnotDb, Description, n.background")
  }
  
  if (!is.numeric(sizecut) || length(sizecut) != 2 || any(is.na(sizecut))) {
    stop("SEPA compressor error: sizecut must be numeric length-2 (min,max)")
  }
  
  if (sizecut[1] > sizecut[2]) {
    stop("SEPA compressor error: sizecut[1] should be <= sizecut[2]")
  }
  
  background_sizes = pa_matrices$terms$n.background
  size_ok = (background_sizes >= sizecut[1]) & (background_sizes <= sizecut[2])
  
  # ============================================================================
  # Combine filters
  # ============================================================================
  
  row_mask = keep_rows & size_ok
  if (!any(row_mask)) {
    stop("SEPA compressor error: no terms remain after filters")
  }
  
  # ============================================================================
  # Apply filters
  # ============================================================================
  
  p_matrix = if (useFDR) pa_matrices$p.adj else pa_matrices$p.raw
  p_subset = p_matrix[row_mask, depth_columns, drop = FALSE]
  n_subset = pa_matrices$n.hits[row_mask, depth_columns, drop = FALSE]
  fc_subset = pa_matrices$fc[row_mask, depth_columns, drop = FALSE]
  
  passed = (p_subset <= pcut) & (n_subset >= ncut) & (fc_subset >= fccut)
  passed[is.na(passed)] = FALSE
  
  if (!any(passed)) {
    return(data.frame(
      term.id = character(0), AnnotDb = character(0), Description = character(0),
      n.background = integer(0),
      is_first = logical(0), is_best = logical(0),
      topn = character(0), p.raw = numeric(0), p.adj = numeric(0), n.hit = integer(0), fc = numeric(0),
      region_start_topn = character(0), region_end_topn = character(0), region_len = integer(0),
      hits.SYMBOL = character(0), hits.ENTREZID = character(0),
      stringsAsFactors = FALSE
    ))
  }
  
  # ============================================================================
  # Prepare annotations
  # ============================================================================
  
  annotations_by_depth = lapply(seq_len(n_depths), function(k) {
    depth_df = do.call(rbind, list(PA$GO[[k]], PA$Reactome[[k]], PA$KEGG[[k]]))
    if (is.null(depth_df) || !nrow(depth_df)) {
      data.frame(term.id = character(0), SYMBOL = character(0), ENTREZID = character(0), stringsAsFactors = FALSE)
    } else {
      required_cols = c("term.id", "SYMBOL", "ENTREZID")
      missing_cols = setdiff(required_cols, names(depth_df))
      for (col in missing_cols) depth_df[[col]] = NA_character_
      depth_df[, required_cols]
    }
  })
  
  # ============================================================================
  # Loop on terms to build output rows
  # ============================================================================
  
  output_rows = vector("list", length = sum(row_mask) * 2L)
  output_index = 0L
  
  kept_term_indices = which(row_mask)
  kept_terms_df = pa_matrices$terms[kept_term_indices, , drop = FALSE]
  kept_term_ids = kept_terms_df$term.id
  
  for (row_index in seq_along(kept_term_ids)) {
    current_term = kept_term_ids[row_index]
    row_passed = passed[row_index, ]
    if (!any(row_passed)) next
    
    regions = find_regions(row_passed, min_len = scut)
    if (!length(regions$starts)) next
    
    # FIRST: start of earliest good region
    rel_first_idx = regions$starts[1]
    abs_first_idx = depth_columns[rel_first_idx]
    
    # BEST: most significant depth
    rel_positions = unlist(mapply(seq.int, regions$starts, regions$ends), use.names = FALSE)
    scores = p_subset[row_index, rel_positions]
    scores = as.numeric(scores)
    scores[!is.finite(scores)] = Inf
    rel_best_idx = rel_positions[which.min(scores)]
    abs_best_idx = depth_columns[rel_best_idx]
    
    term_metadata = kept_terms_df[row_index, , drop = FALSE]
    
    if (abs_best_idx == abs_first_idx) {
      output_index = output_index + 1L
      output_rows[[output_index]] = make_one_row(rel_first_idx, abs_first_idx, TRUE, TRUE)
    } else {
      output_index = output_index + 1L
      output_rows[[output_index]] = make_one_row(rel_first_idx, abs_first_idx, TRUE, FALSE)
      output_index = output_index + 1L
      output_rows[[output_index]] = make_one_row(rel_best_idx, abs_best_idx, FALSE, TRUE)
    }
  }
  
  if (output_index == 0L) {
    return(data.frame(
      term.id = character(0), AnnotDb = character(0), Description = character(0),
      n.background = integer(0),
      is_first = logical(0), is_best = logical(0),
      topn = character(0), p.raw = numeric(0), p.adj = numeric(0), n.hit = integer(0), fc = numeric(0),
      region_start_topn = character(0), region_end_topn = character(0), region_len = integer(0),
      hits.SYMBOL = character(0), hits.ENTREZID = character(0),
      stringsAsFactors = FALSE
    ))
  }
  
  output_df = do.call(rbind, output_rows[seq_len(output_index)])
  rownames(output_df) = NULL
  
  # Order for report: first rows first, then earlier topn, then significance
  position = match(output_df$topn, column_labels)
  order_metric = if (useFDR) output_df$p.adj else output_df$p.raw
  order_indices = order(!output_df$is_first, position, order_metric, output_df$term.id, na.last = TRUE)
  output_df[order_indices, , drop = FALSE]
}
