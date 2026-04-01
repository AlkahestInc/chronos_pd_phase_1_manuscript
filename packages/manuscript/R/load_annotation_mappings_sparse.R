## SPDX-FileCopyrightText: 2026 Alkahest, Inc. <chronospd_data_access@alkahest.com>
##
## SPDX-License-Identifier: MIT

# ------------------------------------------------------------------------------
# Function: load_annotation_mappings_sparse
# Description: Build sparse matrix mapping IDs to pathway terms (GO, KEGG, Reactome)
# Parameters:
#   restrict_ids - optional vector of IDs to restrict to (default: NULL = all)
# Returns: list with sparse matrix M, score_map, and ID mappings
# ------------------------------------------------------------------------------
load_annotation_mappings_sparse <- function(restrict_ids = NULL) {
  # SYMBOL/UNIPROT/ENTREZID --> ENTREZID map
  entrez_ids <- AnnotationDbi::keys(org.Hs.eg.db::org.Hs.eg.db, keytype = "ENTREZID")
  id_map <- data.frame(
    ID_TYPE = "ENTREZID",
    ANNOT = entrez_ids,
    ENTREZID = entrez_ids,
    stringsAsFactors = FALSE
  )
  
  for (id_type in c("SYMBOL", "UNIPROT")) {
    selection <- AnnotationDbi::select(
      org.Hs.eg.db::org.Hs.eg.db,
      keys = AnnotationDbi::keys(org.Hs.eg.db::org.Hs.eg.db, keytype = id_type),
      keytype = id_type,
      columns = "ENTREZID"
    )
    selection <- dplyr::filter(selection, !is.na(ENTREZID))
    id_map <- dplyr::bind_rows(
      id_map,
      data.frame(
        ID_TYPE = id_type,
        ANNOT = selection[[id_type]],
        ENTREZID = selection$ENTREZID,
        stringsAsFactors = FALSE
      )
    )
  }
  
  if (!is.null(restrict_ids)) {
    id_map <- dplyr::filter(id_map, ANNOT %in% restrict_ids)
  }
  
  # Pathways: GO (direct mappings only, no ancestor propagation)
  go_table <- AnnotationDbi::toTable(org.Hs.eg.db::org.Hs.egGO)
  go_table <- dplyr::filter(go_table, gene_id %in% id_map$ENTREZID)
  go_table <- dplyr::mutate(go_table, TERM = paste0("GO:", go_id))
  go_df <- dplyr::select(go_table, ENTREZID = gene_id, TERM)
  
  # Pathways: KEGG
  kegg_links <- KEGGREST::keggLink("pathway", "hsa")
  kegg_df <- data.frame(
    ENTREZID = sub("^hsa:", "", names(kegg_links)),
    TERM = paste0("KEGG:", sub("^path:", "", kegg_links)),
    stringsAsFactors = FALSE
  )
  kegg_df <- dplyr::filter(kegg_df, ENTREZID %in% id_map$ENTREZID)
  
  # Pathways: Reactome
  reactome_to_entrez <- as.list(reactome.db::reactomePATHID2EXTID)
  reactome_df <- do.call(rbind, lapply(names(reactome_to_entrez), function(pathway_id) {
    data.frame(
      ENTREZID = reactome_to_entrez[[pathway_id]],
      TERM = paste0("RX:", pathway_id),
      stringsAsFactors = FALSE
    )
  }))
  reactome_df <- dplyr::filter(reactome_df, ENTREZID %in% id_map$ENTREZID)
  
  # Stack annotation databases, build sparse master term-id matrix
  id_term <- dplyr::bind_rows(
    dplyr::inner_join(dplyr::select(id_map, ANNOT, ENTREZID), go_df, by = "ENTREZID"),
    dplyr::inner_join(dplyr::select(id_map, ANNOT, ENTREZID), kegg_df, by = "ENTREZID"),
    dplyr::inner_join(dplyr::select(id_map, ANNOT, ENTREZID), reactome_df, by = "ENTREZID")
  )
  id_term <- dplyr::distinct(id_term, ANNOT, TERM)
  
  rows <- factor(id_term$ANNOT)
  cols <- factor(id_term$TERM)
  sparse_matrix <- Matrix::sparseMatrix(
    i = as.integer(rows),
    j = as.integer(cols),
    x = 1,
    dimnames = list(levels(rows), levels(cols))
  )
  score_map <- Matrix::rowSums(sparse_matrix)
  
  # Bidirectional mappings
  orig_to_entrez <- split(id_map$ENTREZID, id_map$ANNOT)
  entrez_to_orig <- split(id_map$ANNOT, id_map$ENTREZID)
  
  # Back-maps for runPAs_parallel()
  # SYMBOL
  symbol_selection <- AnnotationDbi::select(
    org.Hs.eg.db::org.Hs.eg.db,
    keys = unique(id_map$ENTREZID),
    keytype = "ENTREZID",
    columns = "SYMBOL"
  )
  symbol_map <- split(symbol_selection$SYMBOL, symbol_selection$ENTREZID)
  
  # UNIPROT
  uniprot_selection <- AnnotationDbi::select(
    org.Hs.eg.db::org.Hs.eg.db,
    keys = unique(id_map$ENTREZID),
    keytype = "ENTREZID",
    columns = "UNIPROT"
  )
  uniprot_map <- split(uniprot_selection$UNIPROT, uniprot_selection$ENTREZID)
  
  # Output
  list(
    M = sparse_matrix,
    score_map = score_map,
    orig2ent = orig_to_entrez,
    entrez2orig = entrez_to_orig,
    SYMBOL_map = symbol_map,
    UNIPROT_map = uniprot_map
  )
}
