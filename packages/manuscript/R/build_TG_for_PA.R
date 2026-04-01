## SPDX-FileCopyrightText: 2026 Alkahest, Inc. <chronospd_data_access@alkahest.com>
##
## SPDX-License-Identifier: MIT

# ------------------------------------------------------------------------------
# Function: build_TG_for_PA
# Description: Build term-gene association matrix mimicking PA package behavior
#              (propagated GO, KEGG/Reactome naming conventions)
# Parameters:
#   restrict_ids - optional vector of IDs to restrict to (default: NULL = all)
#   include_id_types - ID types to include (default: c("ENTREZID","SYMBOL","UNIPROT"))
# Returns: list with sparse matrix M and ID mappings
# ------------------------------------------------------------------------------
build_TG_for_PA <- function(restrict_ids = NULL, include_id_types = c("ENTREZID", "SYMBOL", "UNIPROT")) {
  stopifnot(is.null(restrict_ids) || is.vector(restrict_ids))
  include_id_types = unique(match.arg(include_id_types, c("ENTREZID", "SYMBOL", "UNIPROT"), several.ok = TRUE))
  
  if (!requireNamespace("AnnotationDbi", quietly = TRUE) ||
      !requireNamespace("org.Hs.eg.db", quietly = TRUE) ||
      !requireNamespace("Matrix", quietly = TRUE)) {
    stop("TGmat: need packages AnnotationDbi, org.Hs.eg.db, Matrix")
  }
  
  # ============================================================================
  # Internal helper functions
  # ============================================================================
  
  # GO naming convention
  normalize_go = function(x) {
    x = trimws(as.character(x))
    x = ifelse(grepl("^GO:", x), sub("^GO:", "", x), x)
    paste0("GO:", x)
  }
  
  # Character trimmer
  trim_char = function(x) trimws(as.character(x))
  
  # ============================================================================
  # Build ANNOT --> ENTREZID universe
  # ============================================================================
  
  entrez_ids = AnnotationDbi::keys(org.Hs.eg.db::org.Hs.eg.db, keytype = "ENTREZID")
  id_map = data.frame(
    ID_TYPE = "ENTREZID",
    ANNOT = entrez_ids,
    ENTREZID = entrez_ids,
    stringsAsFactors = FALSE
  )
  
  # Add other ID types
  for (id_type in intersect(c("SYMBOL", "UNIPROT"), include_id_types)) {
    selection = AnnotationDbi::select(
      org.Hs.eg.db::org.Hs.eg.db,
      keys = AnnotationDbi::keys(org.Hs.eg.db::org.Hs.eg.db, keytype = id_type),
      keytype = id_type,
      columns = "ENTREZID"
    )
    selection = selection[!is.na(selection$ENTREZID), , drop = FALSE]
    id_map = rbind(
      id_map,
      data.frame(
        ID_TYPE = id_type,
        ANNOT = selection[[id_type]],
        ENTREZID = selection$ENTREZID,
        stringsAsFactors = FALSE
      )
    )
  }
  
  # Keep only requested ID types
  id_map = id_map[id_map$ID_TYPE %in% include_id_types, , drop = FALSE]
  
  # Optional restriction to predefined set
  if (!is.null(restrict_ids)) {
    restrict_ids = unique(trim_char(restrict_ids))
    id_map = id_map[id_map$ANNOT %in% restrict_ids, , drop = FALSE]
    if (!nrow(id_map)) stop("TGmat: after restricting to requested IDs, nothing remains!")
  }
  
  # Clean up
  id_map$ANNOT = trim_char(id_map$ANNOT)
  id_map$ENTREZID = trim_char(id_map$ENTREZID)
  keep_entrez = unique(id_map$ENTREZID)
  
  # ============================================================================
  # Build TERM -> ENTREZID mappings (mimic PA package logic)
  # ============================================================================
  
  # GO (propagated with ancestor propagation)
  go_to_all = AnnotationDbi::toTable(org.Hs.eg.db::org.Hs.egGO2ALLEGS)
  go_column = intersect(c("go_id", "GOID"), colnames(go_to_all))[1]
  gene_column = intersect(c("gene_id", "GeneID"), colnames(go_to_all))[1]
  if (is.na(go_column) || is.na(gene_column)) {
    stop("TGmat: unexpected columns in org.Hs.egGO2ALLEGS")
  }
  go_to_all[[go_column]] = normalize_go(go_to_all[[go_column]])
  go_to_all[[gene_column]] = trim_char(go_to_all[[gene_column]])
  go_to_all = go_to_all[go_to_all[[gene_column]] %in% keep_entrez, , drop = FALSE]
  go_df = data.frame(
    ENTREZID = go_to_all[[gene_column]],
    TERM = go_to_all[[go_column]],
    stringsAsFactors = FALSE
  )
  
  # KEGG (hsa == human only)
  kegg_df = data.frame(ENTREZID = character(0), TERM = character(0), stringsAsFactors = FALSE)
  if (requireNamespace("KEGGREST", quietly = TRUE)) {
    kegg_links = KEGGREST::keggLink("pathway", "hsa")
    if (length(kegg_links)) {
      entrez_kegg = sub("^hsa:", "", names(kegg_links))
      pathway_kegg = sub("^path:", "", as.character(kegg_links))
      pathway_kegg = sub("^map", "hsa", pathway_kegg)
      keep_kegg = grepl("^hsa\\d+$", pathway_kegg)
      kegg_df = data.frame(
        ENTREZID = trim_char(entrez_kegg[keep_kegg]),
        TERM = trim_char(pathway_kegg[keep_kegg]),
        stringsAsFactors = FALSE
      )
      kegg_df = kegg_df[kegg_df$ENTREZID %in% keep_entrez, , drop = FALSE]
    }
  }
  
  # Reactome (R-HSA-XXXXX)
  reactome_df = data.frame(ENTREZID = character(0), TERM = character(0), stringsAsFactors = FALSE)
  if (requireNamespace("reactome.db", quietly = TRUE)) {
    reactome_to_entrez = as.list(reactome.db::reactomePATHID2EXTID)
    if (length(reactome_to_entrez)) {
      reactome_df = do.call(rbind, lapply(names(reactome_to_entrez), function(pathway_id) {
        if (is.null(reactome_to_entrez[[pathway_id]]) || !length(reactome_to_entrez[[pathway_id]])) return(NULL)
        data.frame(
          ENTREZID = trim_char(reactome_to_entrez[[pathway_id]]),
          TERM = trim_char(pathway_id),
          stringsAsFactors = FALSE
        )
      }))
      reactome_df = reactome_df[reactome_df$ENTREZID %in% keep_entrez, , drop = FALSE]
    }
  }
  
  # ============================================================================
  # Build ANNOT --[ENTREZID]--> TERM relationships
  # ============================================================================
  
  annot_entrez = id_map[, c("ANNOT", "ENTREZID")]
  go_annot_term = merge(annot_entrez, go_df, by = "ENTREZID", all = FALSE)[, c("ANNOT", "TERM")]
  kegg_annot_term = if (nrow(kegg_df)) {
    merge(annot_entrez, kegg_df, by = "ENTREZID", all = FALSE)[, c("ANNOT", "TERM")]
  } else {
    data.frame(ANNOT = character(0), TERM = character(0))
  }
  reactome_annot_term = if (nrow(reactome_df)) {
    merge(annot_entrez, reactome_df, by = "ENTREZID", all = FALSE)[, c("ANNOT", "TERM")]
  } else {
    data.frame(ANNOT = character(0), TERM = character(0))
  }
  
  id_term = unique(rbind(go_annot_term, kegg_annot_term, reactome_annot_term))
  
  if (!nrow(id_term)) {
    stop("TGmat: No ANNOT --> TERM found; check packages, db obsolescence, ID overlap")
  }
  
  # ============================================================================
  # Build master sparse matrix (rows = ANNOT, cols = TERM)
  # ============================================================================
  
  rows = factor(id_term$ANNOT)
  cols = factor(id_term$TERM)
  sparse_matrix = Matrix::sparseMatrix(
    i = as.integer(rows),
    j = as.integer(cols),
    x = 1L,
    dimnames = list(levels(rows), levels(cols))
  )
  
  # ============================================================================
  # Build auxiliary maps
  # ============================================================================
  
  score_map = Matrix::rowSums(sparse_matrix)
  orig_to_entrez = split(id_map$ENTREZID, id_map$ANNOT)
  entrez_to_orig = split(id_map$ANNOT, id_map$ENTREZID)
  
  # ENTREZID to SYMBOL
  symbol_selection = AnnotationDbi::select(
    org.Hs.eg.db::org.Hs.eg.db,
    keys = unique(id_map$ENTREZID),
    keytype = "ENTREZID",
    columns = "SYMBOL"
  )
  symbol_map = split(symbol_selection$SYMBOL, symbol_selection$ENTREZID)
  
  # ENTREZID to UNIPROT
  uniprot_selection = AnnotationDbi::select(
    org.Hs.eg.db::org.Hs.eg.db,
    keys = unique(id_map$ENTREZID),
    keytype = "ENTREZID",
    columns = "UNIPROT"
  )
  uniprot_map = split(uniprot_selection$UNIPROT, uniprot_selection$ENTREZID)
  
  # ============================================================================
  # Output
  # ============================================================================
  
  list(
    M = sparse_matrix,
    score_map = score_map,
    orig2ent = orig_to_entrez,
    entrez2orig = entrez_to_orig,
    SYMBOL_map = symbol_map,
    UNIPROT_map = uniprot_map
  )
}
