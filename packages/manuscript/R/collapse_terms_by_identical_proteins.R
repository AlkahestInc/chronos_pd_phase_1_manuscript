## SPDX-FileCopyrightText: 2026 Alkahest, Inc. <chronospd_data_access@alkahest.com>
##
## SPDX-License-Identifier: MIT

# ------------------------------------------------------------------------------
# Function: collapse_terms_by_identical_proteins
# Description: Network compressor - creates meta-nodes from terms by grouping term sets
#              with same protein connections. Can handle different node attributes separately.
# Parameters:
#   graph - igraph object
#   type_attr - vertex attribute defining "term" or "protein" class (default: "type")
#   term_value - flag for term class (default: "term")
#   protein_value - flag for protein class (default: "protein")
#   label_attr - vertex label attribute (default: "label")
#   keep_members_attr - attribute to store merged members (default: "merged_members")
#   representative - how to name metanodes (default: "least_words")
#   concat_sep - separator for "concat" mode (default: " | ")
# Returns: compressed igraph object
# ------------------------------------------------------------------------------
collapse_terms_by_identical_proteins <- function(
  graph,
  type_attr = "type",
  term_value = "term",
  protein_value = "protein",
  label_attr = "label",
  keep_members_attr = "merged_members",
  representative = c("least_words", "most_words", "shortest_chars", "longest_chars", "concat"),
  concat_sep = " | "
) {
  
  # Basic checks
  stopifnot(inherits(graph, "igraph"))
  if (!(type_attr %in% igraph::vertex_attr_names(graph))) {
    stop(sprintf("graph collapse: vertex attribute '%s' not found.", type_attr))
  }
  
  # ============================================================================
  # Internal helper functions
  # ============================================================================
  
  # Stringify all attributes for equality check
  attribute_signature = function(vertex_idx) {
    values = lapply(vertex_attribute_names, function(attr_name) {
      value = igraph::vertex_attr(graph, attr_name, index = vertex_idx)
      if (length(value) == 0 || all(is.na(value))) {
        "<NA>"
      } else if (length(value) == 1) {
        as.character(ifelse(is.na(value), "<NA>", value))
      } else {
        paste(as.character(replace(value, is.na(value), "<NA>")), collapse = "||")
      }
    })
    paste(paste0(vertex_attribute_names, "=", values), collapse = ";;")
  }
  
  # Sorted protein-neighbor signature
  protein_signature = function(vertex_idx) {
    neighbors = igraph::neighbors(graph, vertex_idx, mode = "all")
    if (!length(neighbors)) return("")
    neighbor_names = vertices$name[as.integer(neighbors)]
    paste(sort(intersect(neighbor_names, protein_vertex_names)), collapse = "|")
  }
  
  # Word count helper
  word_count = function(text) {
    if (is.na(text) || !nzchar(text)) return(Inf)
    length(strsplit(trimws(text), "\\s+")[[1]])
  }
  
  # Choose representative label
  choose_representative = function(row_indices) {
    names_in_group = term_df$name[row_indices]
    labels_in_group = term_df$label[row_indices]
    
    # Concatenate all
    if (representative == "concat") {
      rep_name = sort(names_in_group, decreasing = FALSE)[1]
      other_names = setdiff(names_in_group, rep_name)
      labels_sorted = labels_in_group[order(labels_in_group, names_in_group)]
      concat_label = paste(unique(labels_sorted), collapse = concat_sep)
      return(list(rep = rep_name, others = other_names, concat_label = concat_label))
    }
    
    # Choose one based on rules
    word_counts = vapply(labels_in_group, word_count, numeric(1))
    char_counts = nchar(labels_in_group, type = "chars", allowNA = TRUE, keepNA = TRUE)
    char_counts[is.na(char_counts)] = Inf
    
    order_indices = switch(
      representative,
      least_words = order(word_counts, char_counts, names_in_group),
      most_words = order(-word_counts, char_counts, names_in_group),
      shortest_chars = order(char_counts, names_in_group),
      longest_chars = order(-char_counts, names_in_group)
    )
    
    rep_name = names_in_group[order_indices][1]
    other_names = setdiff(names_in_group, rep_name)
    list(rep = rep_name, others = other_names, concat_label = NULL)
  }
  
  representative = match.arg(representative)
  
  # ============================================================================
  # Get vertices
  # ============================================================================
  
  vertices = igraph::V(graph)
  vertex_names = vertices$name
  vertex_types = igraph::vertex_attr(graph, type_attr)
  term_indices = which(vertex_types == term_value)
  protein_indices = which(vertex_types == protein_value)
  
  if (!length(term_indices) || !length(protein_indices)) {
    warning("graph collapse: no term or no protein vertices; nothing to compress, returning original graph.")
    return(graph)
  }
  
  protein_vertex_names = vertex_names[protein_indices]
  vertex_attribute_names = setdiff(igraph::vertex_attr_names(graph), c("name", label_attr))
  
  # ============================================================================
  # Build term table
  # ============================================================================
  
  term_rows = lapply(term_indices, function(i) {
    label_value = igraph::vertex_attr(graph, label_attr, index = i)
    if (is.null(label_value) || is.na(label_value) || !nzchar(label_value)) {
      label_value = vertex_names[i]
    }
    list(
      idx = i,
      name = vertex_names[i],
      label = as.character(label_value),
      prot_sig = protein_signature(i),
      attr_sig = attribute_signature(i)
    )
  })
  
  term_df = do.call(rbind, lapply(term_rows, as.data.frame, stringsAsFactors = FALSE))
  
  # ============================================================================
  # Group terms by signature
  # ============================================================================
  
  signature_key = paste(term_df$prot_sig, term_df$attr_sig, sep = "###")
  groups = split(seq_len(nrow(term_df)), signature_key)
  
  # Work only with groups size > 1 and non-empty signature
  mergeable_keys = names(groups)[vapply(groups, length, 1L) > 1 & names(groups) != "###"]
  
  if (!length(mergeable_keys)) {
    if ("shape" %in% igraph::vertex_attr_names(graph)) {
      shapes = igraph::vertex_attr(graph, "shape")
      shapes[is.na(shapes) | shapes == ""] = "circle"
      graph = igraph::set_vertex_attr(graph, "shape", value = shapes)
    }
    return(graph)
  }
  
  # ============================================================================
  # Identify vertices to delete and build member map
  # ============================================================================
  
  vertices_to_delete = character(0)
  members_map = list()
  label_overrides = list()
  
  for (key in mergeable_keys) {
    choice = choose_representative(groups[[key]])
    if (length(choice$others)) {
      vertices_to_delete = c(vertices_to_delete, choice$others)
      if (!is.null(keep_members_attr)) {
        members_map[[choice$rep]] = c(members_map[[choice$rep]], choice$others)
      }
      if (!is.null(choice$concat_label)) {
        label_overrides[[choice$rep]] = choice$concat_label
      }
    }
  }
  
  vertices_to_delete = unique(vertices_to_delete)
  if (!length(vertices_to_delete)) {
    if ("shape" %in% igraph::vertex_attr_names(graph)) {
      shapes = igraph::vertex_attr(graph, "shape")
      shapes[is.na(shapes) | shapes == ""] = "circle"
      graph = igraph::set_vertex_attr(graph, "shape", value = shapes)
    }
    return(graph)
  }
  
  # ============================================================================
  # Attach merged members to representative nodes
  # ============================================================================
  
  if (!is.null(keep_members_attr) && length(members_map)) {
    for (rep_name in names(members_map)) {
      new_members = paste(unique(members_map[[rep_name]]), collapse = "|")
      current_value = igraph::vertex_attr(graph, keep_members_attr, index = rep_name)
      if (is.null(current_value) || is.na(current_value) || !nzchar(current_value)) {
        new_value = new_members
      } else {
        new_value = paste(current_value, new_members, sep = "|")
      }
      graph = igraph::set_vertex_attr(graph, keep_members_attr, index = rep_name, value = new_value)
    }
  }
  
  # ============================================================================
  # Override labels for concat mode
  # ============================================================================
  
  if (length(label_overrides)) {
    for (rep_name in names(label_overrides)) {
      graph = igraph::set_vertex_attr(graph, label_attr, index = rep_name, value = label_overrides[[rep_name]])
    }
  }
  
  # ============================================================================
  # Delete non-representatives and simplify
  # ============================================================================
  
  compressed_graph = igraph::delete_vertices(graph, vertices_to_delete)
  compressed_graph = igraph::simplify(compressed_graph, remove.multiple = TRUE, remove.loops = TRUE)
  
  # ============================================================================
  # Vertex shape sanity
  # ============================================================================
  
  if ("shape" %in% igraph::vertex_attr_names(compressed_graph)) {
    shapes = igraph::vertex_attr(compressed_graph, "shape")
    shapes[is.na(shapes) | shapes == ""] = "circle"
    compressed_graph = igraph::set_vertex_attr(compressed_graph, "shape", value = shapes)
  }
  
  compressed_graph
}
