## SPDX-FileCopyrightText: 2026 Alkahest, Inc. <chronospd_data_access@alkahest.com>
##
## SPDX-License-Identifier: MIT

# ------------------------------------------------------------------------------
# Function: runPAs_parallel
# Description: Run GO, KEGG, and Reactome pathway enrichment analysis in parallel
# Parameters:
#   entities - list of character vectors with gene/protein IDs (in id_type format)
#   universe - character vector with background IDs (platform/background in id_type format)
#   id_type - one of "SYMBOL", "ENTREZID", "UNIPROT"
#   BPPARAM - BiocParallel backend configuration
#   ont_go - GO ontology: "BP", "MF", "CC", or "ALL"
#   pval_cutoff - p-value cutoff for enrichment results
#   qval_cutoff - adjusted p-value cutoff for enrichment results
#   minGSSize - minimum gene set size for testing
#   maxGSSize - maximum gene set size for testing
# Returns: list with GO, Reactome, KEGG enrichment results and mapping info
# ------------------------------------------------------------------------------
runPAs_parallel <- function(
    entities,
    universe,
    id_type = c("SYMBOL", "ENTREZID", "UNIPROT"),
    BPPARAM = BiocParallel::SnowParam(workers = max(1, parallel::detectCores() - 1)),
    ont_go = "ALL",
    pval_cutoff = 1,
    qval_cutoff = 1,
    minGSSize = 1,
    maxGSSize = 1e9
) {
  # Required by BiocParallel::bplapply for workers
  library(org.Hs.eg.db)
  library(clusterProfiler)
  library(ReactomePA)
      
  id_type <- match.arg(id_type)
  
  # --- ID mapping: lean + vectorized ---
  if (id_type == "ENTREZID") {
    map_df <- data.frame(
      ENTREZID = universe,
      ORIGINAL = universe,
      stringsAsFactors = FALSE
    )
  } else {
    map_df <- AnnotationDbi::select(
      AnnotationDbi::get("org.Hs.eg.db"),
      keys = unique(c(universe, unlist(entities))),
      keytype = id_type,
      columns = "ENTREZID"
    )
  }
  
  # Vectorized lookup (1st ENTREZID hit or NA)
  map_lookup <- split(map_df$ENTREZID, map_df[[id_type]])
  map_one <- function(x) {
    v <- map_lookup[[x]]
    if (length(v)) v[[1]] else NA_character_
  }
  map_vec <- function(v) {
    vapply(v, map_one, FUN.VALUE = character(1L), USE.NAMES = FALSE)
  }
  
  # Mapping tables
  universe_df <- data.frame(
    entity = universe,
    ENTREZID = map_vec(universe),
    stringsAsFactors = FALSE
  )
  
  entities_df_list <- lapply(entities, function(gene_set) {
    data.frame(
      entity = gene_set,
      ENTREZID = map_vec(gene_set),
      stringsAsFactors = FALSE
    )
  })
  
  # Build mappings for enrichment
  orig_to_entrez <- split(map_df$ENTREZID, map_df[[id_type]])
  symbol_map_df <- AnnotationDbi::select(
    AnnotationDbi::get("org.Hs.eg.db"),
    unique(map_df$ENTREZID),
    "SYMBOL",
    "ENTREZID"
  )
  symbol_map <- split(symbol_map_df$SYMBOL, symbol_map_df$ENTREZID)
  
  # Pre-flatten SYMBOLs per ENTREZID to a single string
  symbol_join <- setNames(
    vapply(
      names(symbol_map),
      function(k) paste(na.omit(symbol_map[[k]]), collapse = "/"),
      FUN.VALUE = character(1L)
    ),
    names(symbol_map)
  )
  
  entrez_sets <- lapply(entities, function(v) {
    unique(unlist(orig_to_entrez[v], use.names = FALSE))
  })
  background_entrez <- unique(unlist(orig_to_entrez[universe], use.names = FALSE))
  
  # ------------------------------------------------------------------------------
  # Internal helper: make_table
  # Description: Build minimal result table from enrichment output
  # ------------------------------------------------------------------------------
  make_table <- function(result_df, annot_db, background_entrez, input_size) {
    if (is.null(result_df) || !nrow(result_df)) {
      return(data.frame())
    }
    
    bg_count <- if (!is.null(result_df$BgRatio)) {
      as.integer(sapply(strsplit(result_df$BgRatio, "/"), `[`, 1))
    } else if (!is.null(result_df$BgCount)) {
      result_df$BgCount
    } else {
      stop("No BgRatio/BgCount")
    }
    
    if (!is.null(result_df$GeneRatio)) {
      gene_ratio <- strsplit(result_df$GeneRatio, "/")
      hit_count <- as.integer(sapply(gene_ratio, "[", 1))
      input_count <- as.integer(sapply(gene_ratio, "[", 2))
    } else {
      hit_count <- result_df$Count
      input_count <- rep(input_size, length(hit_count))
    }
    
    fold_change <- (hit_count / input_count) / (bg_count / length(background_entrez))
    entrez_lists <- strsplit(result_df$geneID, "/")
    entrez_string <- vapply(
      entrez_lists,
      function(e) paste(na.omit(e), collapse = "/"),
      character(1L)
    )
    
    # Main output (SYMBOL filled later on master)
    data.frame(
      term.id = result_df$ID,
      AnnotDb = annot_db,
      Description = result_df$Description,
      n.background = bg_count,
      n.input = input_count,
      n.hit = hit_count,
      FoldChange = fold_change,
      p.value = result_df$pvalue,
      fdr = result_df$p.adjust,
      SYMBOL = NA_character_,
      ENTREZID = entrez_string,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  }
  
  # --- GO enrichment (parallel) ---
  GO_results <- BiocParallel::bplapply(
    X = entrez_sets,
    FUN = function(gene_set, background_entrez, ont_go, pval_cutoff, qval_cutoff, minGSSize, maxGSSize) {
      # Required by BiocParallel::bplapply for workers
      library(org.Hs.eg.db)
      library(clusterProfiler)

      if (length(gene_set) < 2) {
        return(data.frame())
      }
      
      result_df <- as.data.frame(clusterProfiler::enrichGO(
        gene = gene_set,
        universe = background_entrez,
        OrgDb = "org.Hs.eg.db",
        keyType = "ENTREZID",
        ont = ont_go,
        pvalueCutoff = pval_cutoff,
        qvalueCutoff = qval_cutoff,
        minGSSize = minGSSize,
        maxGSSize = maxGSSize,
        readable = FALSE
      ))
      
      annot_db_name <- if (nrow(result_df)) {
        paste0("GO_", result_df$ONTOLOGY)
      } else {
        paste0("GO_", ont_go)
      }
      
      make_table(result_df, annot_db_name, background_entrez, length(gene_set))
    },
    background_entrez = background_entrez,
    ont_go = ont_go,
    pval_cutoff = pval_cutoff,
    qval_cutoff = qval_cutoff,
    minGSSize = minGSSize,
    maxGSSize = maxGSSize,
    BPPARAM = BPPARAM
  )
  
  # --- Reactome enrichment (parallel) ---
  Reactome_results <- BiocParallel::bplapply(
    X = entrez_sets,
    FUN = function(gene_set, background_entrez, pval_cutoff, qval_cutoff, minGSSize, maxGSSize) {
      # Required by BiocParallel::bplapply for workers
      library(ReactomePA)
      
      if (length(gene_set) < 2) {
        return(data.frame())
      }
      
      result_df <- as.data.frame(ReactomePA::enrichPathway(
        gene = gene_set,
        organism = "human",
        universe = background_entrez,
        pvalueCutoff = pval_cutoff,
        qvalueCutoff = qval_cutoff,
        minGSSize = minGSSize,
        maxGSSize = maxGSSize,
        readable = FALSE
      ))
      
      make_table(result_df, "Reactome", background_entrez, length(gene_set))
    },
    background_entrez = background_entrez,
    pval_cutoff = pval_cutoff,
    qval_cutoff = qval_cutoff,
    minGSSize = minGSSize,
    maxGSSize = maxGSSize,
    BPPARAM = BPPARAM
  )
  
  # --- KEGG enrichment (parallel) ---
  KEGG_results <- BiocParallel::bplapply(
    X = entrez_sets,
    FUN = function(gene_set, background_entrez, pval_cutoff, qval_cutoff, minGSSize, maxGSSize) {
      # Required by BiocParallel::bplapply for workers
      library(clusterProfiler)
      
      if (length(gene_set) < 2) {
        return(data.frame())
      }
      
      result_df <- as.data.frame(clusterProfiler::enrichKEGG(
        gene = gene_set,
        organism = "hsa",
        universe = background_entrez,
        pvalueCutoff = pval_cutoff,
        qvalueCutoff = qval_cutoff,
        minGSSize = minGSSize,
        maxGSSize = maxGSSize
      ))
      
      make_table(result_df, "KEGG", background_entrez, length(gene_set))
    },
    background_entrez = background_entrez,
    pval_cutoff = pval_cutoff,
    qval_cutoff = qval_cutoff,
    minGSSize = minGSSize,
    maxGSSize = maxGSSize,
    BPPARAM = BPPARAM
  )
  
  # --- Add SYMBOL annotations to results ---
  add_symbols <- function(result_df) {
    if (is.null(result_df) || !nrow(result_df)) {
      return(result_df)
    }
    
    result_df$SYMBOL <- vapply(
      strsplit(result_df$ENTREZID, "/"),
      function(id_vector) {
        symbol_pieces <- symbol_join[id_vector]
        symbol_pieces <- symbol_pieces[!is.na(symbol_pieces) & symbol_pieces != ""]
        paste(symbol_pieces, collapse = "/")
      },
      FUN.VALUE = character(1L)
    )
    result_df
  }
  
  GO_results <- lapply(GO_results, add_symbols)
  Reactome_results <- lapply(Reactome_results, add_symbols)
  KEGG_results <- lapply(KEGG_results, add_symbols)
  
  # Output
  list(
    GO = GO_results,
    Reactome = Reactome_results,
    KEGG = KEGG_results,
    tested_entities = entities_df_list,
    dataset_background = universe_df,
    id_type = id_type
  )
}
