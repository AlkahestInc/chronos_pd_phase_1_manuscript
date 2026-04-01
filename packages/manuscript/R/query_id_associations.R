## SPDX-FileCopyrightText: 2026 Alkahest, Inc. <chronospd_data_access@alkahest.com>
##
## SPDX-License-Identifier: MIT

# ------------------------------------------------------------------------------
# Function: query_id_associations
# Description: Query pathway term associations for given IDs
# Parameters:
#   maps - output from load_annotation_mappings_sparse()
#   ids - character vector of IDs to query
#   kegg_chunk_size - batch size for KEGG API queries (default: 10)
# Returns: data frame with term_id, term_name, term_size, [input_id]
# ------------------------------------------------------------------------------
query_id_associations <- function(maps, ids, kegg_chunk_size = 10) {
  stopifnot(!is.null(maps$M))
  
  sparse_matrix <- maps$M
  ids <- unique(as.character(ids))
  
  empty_output <- function(present_ids) {
    output <- data.frame(
      term_id = character(),
      term_name = character(),
      term_size = integer(),
      input_id = character(),
      stringsAsFactors = FALSE
    )
    attr(output, "not_found") <- setdiff(ids, present_ids)
    output
  }
  
  present_ids <- ids[ids %in% rownames(sparse_matrix)]
  if (!length(present_ids)) {
    return(empty_output(character()))
  }
  
  term_sizes <- Matrix::colSums(sparse_matrix)
  
  # Get terms for each ID
  terms_per_id <- lapply(present_ids, function(id) {
    term_cols <- which(sparse_matrix[id, ] != 0)
    if (!length(term_cols)) character() else colnames(sparse_matrix)[term_cols]
  })
  names(terms_per_id) <- present_ids
  
  # Convert to long format
  long_format <- do.call(rbind, lapply(names(terms_per_id), function(id) {
    if (!length(terms_per_id[[id]])) {
      return(NULL)
    }
    data.frame(input_id = id, term_id = terms_per_id[[id]], stringsAsFactors = FALSE)
  }))
  if (is.null(long_format) || !nrow(long_format)) {
    return(empty_output(present_ids))
  }
  
  unique_terms <- unique(long_format$term_id)
  is_go <- grepl("^GO:", unique_terms)
  is_kegg <- grepl("^KEGG:", unique_terms)
  is_reactome <- grepl("^RX:", unique_terms)
  
  # --- GO names ---
  go_terms <- unique_terms[is_go]
  go_name_map <- data.frame(
    term_id = go_terms,
    term_name = rep(NA_character_, length(go_terms)),
    stringsAsFactors = FALSE
  )
  if (length(go_terms) && requireNamespace("GO.db", quietly = TRUE)) {
    go_keys <- sub("^GO:", "", go_terms)
    go_selection <- try(
      suppressMessages(
        AnnotationDbi::select(GO.db::GO.db, keys = go_keys, keytype = "GOID", columns = "TERM")
      ),
      silent = TRUE
    )
    if (!inherits(go_selection, "try-error") && nrow(go_selection)) {
      go_selection <- unique(go_selection[, c("GOID", "TERM")])
      idx <- match(sub("^GO:", "", go_name_map$term_id), go_selection$GOID)
      go_name_map$term_name <- go_selection$TERM[idx]
    }
  }
  
  # --- Reactome names ---
  reactome_terms <- unique_terms[is_reactome]
  reactome_name_map <- data.frame(
    term_id = reactome_terms,
    term_name = rep(NA_character_, length(reactome_terms)),
    stringsAsFactors = FALSE
  )
  if (length(reactome_terms) && requireNamespace("reactome.db", quietly = TRUE)) {
    reactome_list <- try(as.list(reactome.db::reactomePATHID2NAME), silent = TRUE)
    if (!inherits(reactome_list, "try-error")) {
      reactome_keys <- sub("^RX:", "", reactome_terms)
      reactome_names <- reactome_list[reactome_keys]
      reactome_names <- reactome_names[!vapply(reactome_names, is.null, logical(1))]
      if (length(reactome_names)) {
        temp_df <- data.frame(
          term_id = paste0("RX:", names(reactome_names)),
          term_name = unname(unlist(reactome_names, use.names = FALSE)),
          stringsAsFactors = FALSE
        )
        reactome_name_map$term_name <- temp_df$term_name[match(reactome_name_map$term_id, temp_df$term_id)]
      }
    }
  }
  
  # --- KEGG names ---
  kegg_terms <- unique_terms[is_kegg]
  kegg_name_map <- data.frame(
    term_id = kegg_terms,
    term_name = rep(NA_character_, length(kegg_terms)),
    stringsAsFactors = FALSE
  )
  if (length(kegg_terms) && requireNamespace("KEGGREST", quietly = TRUE)) {
    kegg_ids <- sub("^KEGG:", "", kegg_terms)
    kegg_chunks <- split(kegg_ids, ceiling(seq_along(kegg_ids) / max(1, kegg_chunk_size)))
    kegg_results <- lapply(kegg_chunks, function(chunk) {
      kegg_data <- try(KEGGREST::keggGet(chunk), silent = TRUE)
      if (inherits(kegg_data, "try-error") || length(kegg_data) == 0) {
        return(NULL)
      }
      do.call(rbind, lapply(kegg_data, function(entry) {
        entry_id <- if (!is.null(entry$ENTRY)) sub("\\s.*$", "", entry$ENTRY) else NA_character_
        entry_name <- if (!is.null(entry$NAME) && length(entry$NAME)) entry$NAME[1] else NA_character_
        entry_name <- sub(" - Homo sapiens \\(.*\\)$", "", entry_name)
        data.frame(kid = entry_id, term_name = entry_name, stringsAsFactors = FALSE)
      }))
    })
    kegg_df <- do.call(rbind, Filter(Negate(is.null), kegg_results))
    if (!is.null(kegg_df) && nrow(kegg_df)) {
      kegg_df$term_id <- paste0("KEGG:", kegg_df$kid)
      kegg_name_map$term_name <- kegg_df$term_name[match(kegg_name_map$term_id, kegg_df$term_id)]
    }
  }
  
  # Bind all name maps
  name_map <- rbind(go_name_map, kegg_name_map, reactome_name_map)
  remaining_terms <- setdiff(unique_terms, name_map$term_id)
  if (length(remaining_terms)) {
    name_map <- rbind(
      name_map,
      data.frame(
        term_id = remaining_terms,
        term_name = NA_character_,
        stringsAsFactors = FALSE
      )
    )
  }
  
  # Final output
  output <- merge(
    long_format,
    data.frame(
      term_id = names(term_sizes),
      term_size = as.integer(term_sizes),
      stringsAsFactors = FALSE
    ),
    by = "term_id",
    all.x = TRUE,
    sort = FALSE
  )
  output <- merge(output, name_map, by = "term_id", all.x = TRUE, sort = FALSE)
  output <- output[order(-output$term_size, output$term_id, output$input_id), ]
  rownames(output) <- NULL
  
  if (length(present_ids) == 1) output$input_id <- NULL
  
  attr(output, "not_found") <- setdiff(ids, present_ids)
  output
}
