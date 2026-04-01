## SPDX-FileCopyrightText: 2026 Alkahest, Inc. <chronospd_data_access@alkahest.com>
##
## SPDX-License-Identifier: MIT

here::i_am("src/20_analyses/sepa/PA_heatmap_network.R")

# ==============================================================================
# SEPA Results Interpreter: Heatmap and Term-Gene Network Visualization
# ==============================================================================
# Author: T.Nanasi
# Date: 2025.11
# Description: Load and visualize pathway enrichment results from SEPA analysis
#              Generate heatmaps and interactive network plots

# ==============================================================================
# ENVIRONMENT SETUP
# ==============================================================================
script_start_time = Sys.time()

# ==============================================================================
# REQUIRED PACKAGES
# ------------------------------------------------------------------------------
# This script requires the following packages:
# - Data manipulation: pracma, dplyr, matrixStats, Matrix
# - Annotation databases: AnnotationDbi, org.Hs.eg.db, KEGGREST, reactome.db
# - Visualization: ComplexHeatmap, circlize, grid, igraph, visNetwork, gplots
# - Clustering: simplifyEnrichment
# - Export (optional): openxlsx
# ==============================================================================

# ==============================================================================
# CONTROL FLAGS - set these to control script behavior
# ==============================================================================
EXPORT_EXCEL = TRUE          # export Excel summary tables
EXPORT_HEATMAP_MAIN = TRUE    # export main heatmap PDF
EXPORT_HEATMAP_SUPPL = TRUE  # export supplementary heatmap PDF
EXPORT_NETWORK_STATIC = TRUE  # export static network PDF

# ==============================================================================
# FILEPATHS - modify these as needed
# ==============================================================================
INPUT_DIR_PA <- here::here(glue::glue("{get_config()$checkpoint_dir}/")) 
OUTPUT_FILE_EXCEL_FIRST_BEST = glue::glue("{get_config()$output_dir}/supplementary_tables/table_s6.xlsx")
dir.create(dirname(OUTPUT_FILE_EXCEL_FIRST_BEST), showWarnings = FALSE, recursive = TRUE)
OUTPUT_FILE_EXCEL_HEATMAP_SUMMARY_MAIN = 'heatmap_term_summary_main.xlsx'
dir.create(dirname(OUTPUT_FILE_EXCEL_HEATMAP_SUMMARY_MAIN), showWarnings = FALSE, recursive = TRUE)
OUTPUT_FILE_EXCEL_HEATMAP_SUMMARY_SUPPL = 'heatmap_term_summary_suppl.xlsx'
dir.create(dirname(OUTPUT_FILE_EXCEL_HEATMAP_SUMMARY_SUPPL), showWarnings = FALSE, recursive = TRUE)
OUTPUT_FILE_HEATMAP_MAIN = glue::glue("{get_config()$output_dir}/main_figures/figure_3/figure_3c.pdf")
dir.create(dirname(OUTPUT_FILE_HEATMAP_MAIN), showWarnings = FALSE, recursive = TRUE)
OUTPUT_FILE_HEATMAP_SUPPL = glue::glue("{get_config()$output_dir}/supplementary_figures/figure_s5/figure_s5a.pdf")
dir.create(dirname(OUTPUT_FILE_HEATMAP_SUPPL), showWarnings = FALSE, recursive = TRUE)
OUTPUT_FILE_NETWORK_STATIC = glue::glue("{get_config()$output_dir}/main_figures/figure_3/figure_3d.pdf")
dir.create(dirname(OUTPUT_FILE_NETWORK_STATIC), showWarnings = FALSE, recursive = TRUE)
OUTPUT_FILE_NETWORK_STATIC_LEGEND = glue::glue("{get_config()$output_dir}/main_figures/figure_3/figure_3d-legend.pdf")
dir.create(dirname(OUTPUT_FILE_NETWORK_STATIC_LEGEND), showWarnings = FALSE, recursive = TRUE)

# ==============================================================================
# CONSTANTS - PA condition names from SEPA_classic.R output
# ==============================================================================
PA_CONDITION_ALL_UP = 'lme_SomaOlinkAlamar_all_up'
PA_CONDITION_ALL_DOWN = 'lme_SomaOlinkAlamar_all_down'
PA_CONDITION_PRE_UP = 'lme_SomaOlinkAlamar_pre_up'
PA_CONDITION_PRE_DOWN = 'lme_SomaOlinkAlamar_pre_down'
PA_CONDITION_POST_UP = 'lme_SomaOlinkAlamar_post_up'
PA_CONDITION_POST_DOWN = 'lme_SomaOlinkAlamar_post_down'
PA_CONDITION_PREJM_UP = 'lme_SomaOlinkAlamar_prejm_up'
PA_CONDITION_PREJM_DOWN = 'lme_SomaOlinkAlamar_prejm_down'
PA_CONDITION_PREJM_MIX = 'lme_SomaOlinkAlamar_prejm_mix'

# ==============================================================================
# HEATMAP FILTER SETTINGS
# ==============================================================================
HEATMAP_USE_FDR = TRUE
HEATMAP_P_CUTOFF = 0.01
HEATMAP_N_CUTOFF = 3
HEATMAP_FC_CUTOFF = 2
HEATMAP_SIZE_CUTOFF = c(30, 2000)
HEATMAP_DEPTH_RANGE = 1:15
HEATMAP_CONSECUTIVE_STEPS = 3
HEATMAP_ANNOTATION_DBS = c('KEGG', 'GO_BP', 'Reactome')
HEATMAP_SHOW_ROW_NAMES = TRUE
HEATMAP_ROW_TITLE_ROT = 0
HEATMAP_RECT_BORDER_COLOR = 'dimgray'
HEATMAP_PRESELECTED_TERMS = NULL

HEATMAP_SELECTED_CONDITIONS = c(
  PA_CONDITION_ALL_UP,
  PA_CONDITION_ALL_DOWN,
  PA_CONDITION_PRE_UP,
  PA_CONDITION_PRE_DOWN,
  PA_CONDITION_POST_UP,
  PA_CONDITION_POST_DOWN
)

# ==============================================================================
# NETWORK PLOT SETTINGS
# ==============================================================================
NETWORK_TOP_N_TERMS = 5   # number of top terms per condition to include in network
NETWORK_GENE_DEPTH = 2    # depth level to extract genes from PA results

NETWORK_CONDITION_NAMES = c(
  PA_CONDITION_PRE_UP,
  PA_CONDITION_PRE_DOWN,
  PA_CONDITION_POST_UP,
  PA_CONDITION_POST_DOWN
)

# ==============================================================================
# PART 1: DATA LOADING & PA RESULTS FORMATTING
# ==============================================================================
{
print('Loading PA results from directory...')

# ------------------------------------------------------------------------------
# Load PA run results from files
# ------------------------------------------------------------------------------

PA_results = list()
pa_files = list.files(
  INPUT_DIR_PA,
  pattern = "SEPA_UNIPROT_.*.Rdata"
)

for (k in 1:length(pa_files)) {
  print(paste('  loading', pa_files[k], '...'))
  temp_env = new.env()
  load(file.path(INPUT_DIR_PA, pa_files[k]), temp_env)
  PA_results[[length(PA_results) + 1]] = temp_env
  names(PA_results)[[length(PA_results)]] = temp_env$used.filters
  rm(temp_env)
}

# ------------------------------------------------------------------------------
# Safety check: are entities encoded the same way?
# ------------------------------------------------------------------------------

mapping_check = compare_uniprot_mappings(PA_results)
print(mapping_check)
if (mean(mapping_check$same_mapping_in_shared) != 1) {
  stop('mapping conflicts in PAs')
}

# ------------------------------------------------------------------------------
# Convert PA results to matrix format
# ------------------------------------------------------------------------------

PA_matrices = list()
for (q in 1:length(PA_results)) {
  PA_matrices[[q]] = summarize_PA_matrix(PA_results[[q]]$PA)
}
names(PA_matrices) = names(PA_results)

# ------------------------------------------------------------------------------
# Sync PA matrices to common term set
# ------------------------------------------------------------------------------

# Combine all terms from all conditions
all_terms_combined = do.call(rbind, lapply(PA_matrices, function(pa_mat) pa_mat$terms))

# Aggregate background sizes (min/max across conditions)
columns_to_group = setdiff(colnames(all_terms_combined), "n.background")
background_min = aggregate(
  all_terms_combined$n.background,
  by = all_terms_combined[columns_to_group],
  FUN = min
)
colnames(background_min)[ncol(background_min)] = "n.background.min"

background_max = aggregate(
  all_terms_combined$n.background,
  by = all_terms_combined[columns_to_group],
  FUN = max
)
colnames(background_max)[ncol(background_max)] = "n.background.max"

# Master list of terms
all_terms = merge(background_min, background_max, by = columns_to_group, sort = FALSE)
rm(all_terms_combined, background_min, background_max, columns_to_group)

# Build synchronized matrices
PA_matrices_synced = list()
for (q in 1:length(PA_matrices)) {
  PA_matrices_synced[[q]] = sync_rows_to_term_ids(PA_matrices[[q]], all_terms$term.id)
}
names(PA_matrices_synced) = names(PA_matrices)

print('PA results loaded and synchronized.')

}  # END PART 1


# ==============================================================================
# PART 2: SUPPLEMENTARY TABLE & ANNOTATION PREPARATION
# ==============================================================================
{
print('Creating supplementary table: global summary of SEPA scan results')

selected_conditions = c(
  PA_CONDITION_PRE_UP,
  PA_CONDITION_PRE_DOWN,
  PA_CONDITION_POST_UP,
  PA_CONDITION_POST_DOWN
)

pa_subset = PA_results[which(names(PA_results) %in% selected_conditions)]
supplementary_table = NULL

for (q in 1:length(pa_subset)) {
  table_rows = build_pa_first_best_table(pa_subset[[q]]$PA)
  if (nrow(table_rows) > 0) {
    table_rows$condition = names(pa_subset)[q]
    table_rows = table_rows[, c(1:4, 17, 5:16)]
    supplementary_table = rbind(supplementary_table, table_rows)
  }
}

# Validate uniqueness of term IDs
unique_check = supplementary_table[, c(1, 4)]
unique_check = unique(unique_check)
if (!length(unique(unique_check$term.id)) == nrow(unique_check)) {
  stop('supplementary table error: non-unique term ids')
}
rm(unique_check, table_rows, pa_subset, selected_conditions)

# Convert numeric columns
supplementary_table$topn = as.numeric(supplementary_table$topn)
supplementary_table$region_start_topn = as.numeric(supplementary_table$region_start_topn)
supplementary_table$region_end_topn = as.numeric(supplementary_table$region_end_topn)

# Export if flag is set
if (EXPORT_EXCEL) {
  if (requireNamespace("openxlsx", quietly = TRUE)) {
    openxlsx::write.xlsx(supplementary_table, OUTPUT_FILE_EXCEL_FIRST_BEST)
    print(paste('Exported:', OUTPUT_FILE_EXCEL_FIRST_BEST))
  } else {
    warning("openxlsx package not available; skipping Excel export")
  }
}


# ==============================================================================
# PREPARE ANNOTATIONS TO MIMIC PA LOGIC
# ==============================================================================

print('Preparing annotations to mimic PA logic...')

annotation_maps = build_TG_for_PA()
term_gene_matrix = annotation_maps$M

# Verify term matching
if (mean(all_terms$term.id %in% colnames(term_gene_matrix)) != 1) {
  stop('term association mismatch - propagated GO terms, deprecated versions..?')
}

print('Annotations prepared.')

}  # END PART 2


# ==============================================================================
# PART 3: HEATMAP GENERATION (MAIN & SUPPLEMENTARY)
# ==============================================================================
{
print('Starting heatmap & network plot generation with key terms')

# Generate both main (clustered) and supplementary (unclustered) figures
for (figure_type in c('main', 'supplementary')) {

  print(paste0('========================================'))
  print(paste0('Generating ', figure_type, ' figure...'))
  print(paste0('========================================'))

  # ============================================================================
  # Apply filter settings from configuration (defined at top of script)
  # ============================================================================

  use_fdr = HEATMAP_USE_FDR
  p_cutoff = HEATMAP_P_CUTOFF
  n_cutoff = HEATMAP_N_CUTOFF
  fc_cutoff = HEATMAP_FC_CUTOFF
  size_cutoff = HEATMAP_SIZE_CUTOFF
  depth_range = HEATMAP_DEPTH_RANGE
  consecutive_steps_cutoff = HEATMAP_CONSECUTIVE_STEPS
  annotation_databases = HEATMAP_ANNOTATION_DBS
  show_row_names = HEATMAP_SHOW_ROW_NAMES
  row_title_rotation = HEATMAP_ROW_TITLE_ROT
  rectangle_border_color = HEATMAP_RECT_BORDER_COLOR
  preselected_terms = HEATMAP_PRESELECTED_TERMS
  selected_pa_conditions = HEATMAP_SELECTED_CONDITIONS

  # Set clustering method and output files based on figure type
  if (figure_type == 'main') {
    term_clustering_method = 'louvain'
    output_excel_file = OUTPUT_FILE_EXCEL_HEATMAP_SUMMARY_MAIN
    output_pdf_file = OUTPUT_FILE_HEATMAP_MAIN
    output_pdf_width = 14
    output_pdf_height = 4
    export_heatmap_flag = EXPORT_HEATMAP_MAIN
  } else {
    term_clustering_method = 'none'
    output_excel_file = OUTPUT_FILE_EXCEL_HEATMAP_SUMMARY_SUPPL
    output_pdf_file = OUTPUT_FILE_HEATMAP_SUPPL
    output_pdf_width = 14
    output_pdf_height = 14
    export_heatmap_flag = EXPORT_HEATMAP_SUPPL
  }

  print(paste('Clustering method:', term_clustering_method))
  
  # ============================================================================
  # Generate SEPA summaries for selected conditions
  # ============================================================================
  
  pa_subset = PA_results[which(names(PA_results) %in% selected_pa_conditions)]
  pa_summaries = list()
  
  for (q in 1:length(pa_subset)) {
    pa_run = pa_subset[[q]]$PA
    pa_run$tested_entities = pa_run$tested_entities[depth_range]
    pa_run$GO = pa_run$GO[depth_range]
    pa_run$KEGG = pa_run$KEGG[depth_range]
    pa_run$Reactome = pa_run$Reactome[depth_range]
    
    summary_df = summarize_significant_terms(
      pa_run,
      pcut = p_cutoff,
      ncut = n_cutoff,
      fccut = fc_cutoff,
      scut = consecutive_steps_cutoff,
      sizecut = size_cutoff,
      useFDR = use_fdr
    )
    
    summary_df = summary_df[which(summary_df$AnnotDb %in% annotation_databases), ]
    
    if (nrow(summary_df) > 0) {
      # Scoring: both early appearance and ultimate strength
      summary_df$score = -log10(summary_df$p.raw_at_best) / as.numeric(summary_df$first_depth_pass)
      summary_df$rank = rank(-summary_df$score)
    }
    
    pa_summaries[[q]] = summary_df
  }
  names(pa_summaries) = names(pa_subset)
  
  # Collect selected terms and genes
  selected_terms = NULL
  selected_genes = NULL
  for (q in 1:length(pa_summaries)) {
    selected_terms = c(selected_terms, pa_summaries[[q]]$term.id)
    selected_genes = unique(rbind(selected_genes, PA_results[[q]]$PA$tested_entities[[max(depth_range)]]))
  }
  selected_terms = unique(selected_terms)
  selected_terms = all_terms$term.id[which(
    all_terms$term.id %in% selected_terms & all_terms$AnnotDb %in% annotation_databases
  )]
  
  if (length(preselected_terms) != 0) {
    selected_terms = intersect(selected_terms, preselected_terms)
  }
  
  print(paste('Terms to plot:', length(selected_terms)))
  print(paste('Genes in test sets:', nrow(selected_genes)))
  
  # ============================================================================
  # Build term-to-gene sets and determine directionality
  # ============================================================================
  
  term_gene_subset = term_gene_matrix[which(rownames(term_gene_matrix) %in% selected_genes$ENTREZID), ]
  term_gene_subset = term_gene_subset[, which(colnames(term_gene_subset) %in% selected_terms)]
  term_gene_subset = term_gene_subset[which(Matrix::rowSums(term_gene_subset) > 0), ]
  
  print(paste('Active genes:', nrow(term_gene_subset)))
  
  unique_gene_sets = as.data.frame(t(as.matrix(term_gene_subset)))
  print(paste(nrow(unique_gene_sets), 'terms --', nrow(unique(unique_gene_sets)), 'unique gene sets'))
  
  # Directionality catalog
  
  # All samples
  all_up_slot = PA_matrices_synced[[PA_CONDITION_ALL_UP]]
  all_up_indices = which(all_up_slot$terms$term.id %in% selected_terms)
  fdr_all_up = all_up_slot$p.adj[all_up_indices, depth_range]
  rownames(fdr_all_up) = paste(all_terms$AnnotDb[all_up_indices], all_terms$Description[all_up_indices], sep = ' : ')
  
  all_down_slot = PA_matrices_synced[[PA_CONDITION_ALL_DOWN]]
  all_down_indices = which(all_down_slot$terms$term.id %in% selected_terms)
  fdr_all_down = all_down_slot$p.adj[all_down_indices, depth_range]
  rownames(fdr_all_down) = paste(all_terms$AnnotDb[all_down_indices], all_terms$Description[all_down_indices], sep = ' : ')
  
  # Pre diagnosis
  pre_up_slot = PA_matrices_synced[[PA_CONDITION_PRE_UP]]
  pre_up_indices = which(pre_up_slot$terms$term.id %in% selected_terms)
  fdr_pre_up = pre_up_slot$p.adj[pre_up_indices, depth_range]
  rownames(fdr_pre_up) = paste(all_terms$AnnotDb[pre_up_indices], all_terms$Description[pre_up_indices], sep = ' : ')
  
  pre_down_slot = PA_matrices_synced[[PA_CONDITION_PRE_DOWN]]
  pre_down_indices = which(pre_down_slot$terms$term.id %in% selected_terms)
  fdr_pre_down = pre_down_slot$p.adj[pre_down_indices, depth_range]
  rownames(fdr_pre_down) = paste(all_terms$AnnotDb[pre_down_indices], all_terms$Description[pre_down_indices], sep = ' : ')
  
  # Post diagnosis
  post_up_slot = PA_matrices_synced[[PA_CONDITION_POST_UP]]
  post_up_indices = which(post_up_slot$terms$term.id %in% selected_terms)
  fdr_post_up = post_up_slot$p.adj[post_up_indices, depth_range]
  rownames(fdr_post_up) = paste(all_terms$AnnotDb[post_up_indices], all_terms$Description[post_up_indices], sep = ' : ')
  
  post_down_slot = PA_matrices_synced[[PA_CONDITION_POST_DOWN]]
  post_down_indices = which(post_down_slot$terms$term.id %in% selected_terms)
  fdr_post_down = post_down_slot$p.adj[post_down_indices, depth_range]
  rownames(fdr_post_down) = paste(all_terms$AnnotDb[post_down_indices], all_terms$Description[post_down_indices], sep = ' : ')
  
  # Categorize terms by directionality
  up_indices = which(
    matrixStats::rowMins(fdr_all_up, na.rm = TRUE) <= 0.05 |
    matrixStats::rowMins(fdr_pre_up, na.rm = TRUE) <= 0.05 |
    matrixStats::rowMins(fdr_post_up, na.rm = TRUE) <= 0.05
  )
  down_indices = which(
    matrixStats::rowMins(fdr_all_down, na.rm = TRUE) <= 0.05 |
    matrixStats::rowMins(fdr_pre_down, na.rm = TRUE) <= 0.05 |
    matrixStats::rowMins(fdr_post_down, na.rm = TRUE) <= 0.05
  )
  both_indices = intersect(up_indices, down_indices)
  up_indices = setdiff(up_indices, both_indices)
  down_indices = setdiff(down_indices, both_indices)
  
  terms_subset = all_terms[which(all_terms$term.id %in% selected_terms), ]
  selected_terms_up = terms_subset$term.id[up_indices]
  selected_terms_down = terms_subset$term.id[down_indices]
  selected_terms_bidi = terms_subset$term.id[both_indices]
  
  # Clean up temporary variables
  rm(
    all_up_slot, all_down_slot, pre_up_slot, pre_down_slot, post_up_slot, post_down_slot,
    all_up_indices, all_down_indices, pre_up_indices, pre_down_indices, post_up_indices, post_down_indices,
    up_indices, down_indices, both_indices, terms_subset, unique_gene_sets
  )

  
  # ============================================================================
  # Cluster terms
  # ============================================================================
  
  if (identical(term_clustering_method, 'none')) {
    
    # No clustering - each term is its own cluster
    terms_up_df = all_terms[which(all_terms$term.id %in% selected_terms_up), ]
    clustering_up = list()
    clustering_up$clusters = data.frame(term.id = terms_up_df$term.id, cluster = 1:nrow(terms_up_df))
    
    terms_down_df = all_terms[which(all_terms$term.id %in% selected_terms_down), ]
    clustering_down = list()
    clustering_down$clusters = data.frame(term.id = terms_down_df$term.id, cluster = 1:nrow(terms_down_df))
    
    terms_bidi_df = all_terms[which(all_terms$term.id %in% selected_terms_bidi), ]
    clustering_bidi = list()
    clustering_bidi$clusters = data.frame(term.id = terms_bidi_df$term.id, cluster = 1:nrow(terms_bidi_df))
    
    rm(terms_up_df, terms_down_df, terms_bidi_df)
    
  } else {
    
    # Validate clustering method
    if (!term_clustering_method %in% simplifyEnrichment::all_clustering_methods()) {
      stop(paste(
        'selected term.clustering.method [', term_clustering_method,
        '] is not valid. possible options are: [',
        paste(simplifyEnrichment::all_clustering_methods(), collapse = ' / '), ']'
      ))
    }
    
    # Cluster each directionality group
    set.seed(0)
    clustering_up = cluster_terms_simplify(selected_terms_up, term_gene_matrix, method = term_clustering_method)
    set.seed(0)
    clustering_down = cluster_terms_simplify(selected_terms_down, term_gene_matrix, method = term_clustering_method)
    set.seed(0)
    clustering_bidi = cluster_terms_simplify(selected_terms_bidi, term_gene_matrix, method = term_clustering_method)
  }
  
  # Merge clusters to single list
  cluster_up_df = clustering_up$clusters[, 1:2]
  cluster_bidi_df = clustering_bidi$clusters[, 1:2]
  cluster_down_df = clustering_down$clusters[, 1:2]
  
  cluster_bidi_df$cluster = cluster_bidi_df$cluster + max(cluster_up_df$cluster)
  cluster_down_df$cluster = cluster_down_df$cluster + max(cluster_bidi_df$cluster)
  
  combined_clustering = list()
  combined_clustering$clusters = rbind(cluster_up_df, cluster_bidi_df, cluster_down_df)
  
  # Join with term metadata
  term_clusters = dplyr::left_join(combined_clustering$clusters, all_terms, by = 'term.id')
  
  # Track ranks for cluster label selection
  for (q in 1:length(pa_summaries)) {
    if (nrow(pa_summaries[[q]]) > 0) {
      rank_df = data.frame(term.id = pa_summaries[[q]]$term.id, rank = pa_summaries[[q]]$rank)
      term_clusters = dplyr::left_join(term_clusters, rank_df, by = 'term.id')
      colnames(term_clusters)[ncol(term_clusters)] = paste0('rank.', names(pa_summaries)[[q]])
    }
  }
  
  rank_columns = which(colnames(term_clusters) %in% paste0('rank.', names(pa_summaries)))
  term_clusters$best.rank = matrixStats::rowMins(as.matrix(term_clusters[, rank_columns]), na.rm = TRUE)
  
  # Define cluster names (select most representative term per cluster)
  cluster_ids = sort(unique(term_clusters$cluster))
  cluster_names = rep(NA, length(cluster_ids))
  names(cluster_names) = cluster_ids
  
  for (k in 1:length(cluster_names)) {
    cluster_id = as.numeric(names(cluster_names)[k])
    candidate_indices = which(term_clusters$cluster == cluster_id)
    
    # Select by: best rank, smallest background, fewest words, shortest description
    if (length(candidate_indices) > 1) {
      candidate_indices = candidate_indices[which(
        term_clusters$best.rank[candidate_indices] == min(term_clusters$best.rank[candidate_indices])
      )]
    }
    if (length(candidate_indices) > 1) {
      candidate_indices = candidate_indices[which(
        term_clusters$n.background.min[candidate_indices] == min(term_clusters$n.background.min[candidate_indices])
      )]
    }
    if (length(candidate_indices) > 1) {
      word_counts = sapply(strsplit(term_clusters$Description[candidate_indices], ' '), length)
      candidate_indices = candidate_indices[which(word_counts == min(word_counts))]
    }
    if (length(candidate_indices) > 1) {
      char_counts = nchar(term_clusters$Description[candidate_indices])
      candidate_indices = candidate_indices[which(char_counts == min(char_counts))]
    }
    
    if (length(candidate_indices) > 1) {
      stop('cannot select single name for term cluster')
    }
    
    cluster_names[k] = term_clusters$Description[candidate_indices]
  }
  
  # Build supplementary table for visualized terms
  term_summary_df = dplyr::left_join(combined_clustering$clusters, all_terms, by = 'term.id')
  term_summary_df = data.frame(
    term.id = term_summary_df$term.id,
    AnnotDb = term_summary_df$AnnotDb,
    Description = term_summary_df$Description,
    n.background = term_summary_df$n.background.min,
    cluster = term_summary_df$cluster
  )
  
  cluster_label_df = data.frame(cluster = as.numeric(names(cluster_names)), label = cluster_names)
  term_summary_df = dplyr::left_join(term_summary_df, cluster_label_df, by = 'cluster')
  
  term_summary_df$directionality = NA
  term_summary_df$directionality[which(term_summary_df$term.id %in% clustering_up$clusters$term.id)] = 'UP'
  term_summary_df$directionality[which(term_summary_df$term.id %in% clustering_down$clusters$term.id)] = 'DOWN'
  term_summary_df$directionality[which(term_summary_df$term.id %in% clustering_bidi$clusters$term.id)] = 'MIX'
  
  # Add top genes per term
  term_summary_df$top_genes = NA
  for (k in 1:nrow(term_summary_df)) {
    term_id = term_summary_df$term.id[k]
    gene_indices = which(term_gene_subset[, which(colnames(term_gene_subset) == term_id)] > 0)
    if (length(gene_indices) > 0) {
      gene_symbols = map_entrez_to_symbols_uniprot(rownames(term_gene_subset)[gene_indices], annotation_maps)$SYMBOL
      term_summary_df$top_genes[k] = paste(sort(unique(gene_symbols)), collapse = '/')
    }
  }
  
  # Export if flag is set
  if (EXPORT_EXCEL) {
    if (requireNamespace("openxlsx", quietly = TRUE)) {
      openxlsx::write.xlsx(term_summary_df, output_excel_file)
      print(paste('Exported:', output_excel_file))
    } else {
      warning("openxlsx package not available; skipping Excel export")
    }
  }
  
  # Clean up
  rm(cluster_up_df, cluster_bidi_df, cluster_down_df, cluster_label_df, rank_columns, cluster_ids, candidate_indices)
  
  # ============================================================================
  # Draw heatmap
  # ============================================================================

  if (identical(term_clustering_method, 'none')) {
    append_term_count = FALSE
  } else {
    append_term_count = TRUE
  }

  # Export heatmap PDF if flag is set
  if (export_heatmap_flag) {
    pdf(output_pdf_file, width = output_pdf_width, height = output_pdf_height)
    plot_cluster_triple_heatmap(
      sel_terms = selected_terms,
      clustering_result = combined_clustering,
      pa_matrices_synced = PA_matrices_synced,
      drange = depth_range,
      fdr_cut = 0.05,
      cluster_names = cluster_names,
      append_n = append_term_count,
      show_row_names = TRUE,
      rect_gp_col = "grey20",
      row_title_rot = 0
    )
    dev.off()
    print(paste('Exported:', output_pdf_file))
  }

  # ============================================================================
  # Network generation (only for main figure to avoid duplication)
  # ============================================================================

  if (figure_type == 'main') {

    print('========================================')
    print('Generating network plot...')
    print('========================================')

    # Apply network settings from configuration
    top_n_terms = NETWORK_TOP_N_TERMS
    gene_depth = NETWORK_GENE_DEPTH
    network_condition_names = NETWORK_CONDITION_NAMES

    # --------------------------------------------------------------------------
    # Select terms for network plot
    # --------------------------------------------------------------------------

    # Select top N terms per condition
    n_top = top_n_terms
    net_terms = NULL
    for (condition_idx in 1:length(network_condition_names)) {
      condition_name = network_condition_names[condition_idx]
      condition_pasum = pa_summaries[[which(names(pa_summaries) == condition_name)]]
      # Filter to only terms shown on heatmap
      condition_pasum = condition_pasum[which(condition_pasum$Description %in% cluster_names), ]
      condition_pasum$rank = rank(condition_pasum$rank)
      net_terms = unique(c(net_terms, condition_pasum$term.id[which(condition_pasum$rank <= n_top)]))
    }
    rm(condition_idx, condition_name, condition_pasum, n_top)

    # Join with all terms info
    net_terms = dplyr::left_join(data.frame(term.id = net_terms), all_terms, by = 'term.id')

    # Mark which conditions each term appears in
    for (condition_idx in 1:length(network_condition_names)) {
      condition_name = network_condition_names[condition_idx]
      condition_pasum = pa_summaries[[which(names(pa_summaries) == condition_name)]]
      net_terms$condition_flag = 0
      net_terms$condition_flag[which(net_terms$term.id %in% condition_pasum$term.id)] = 1
      colnames(net_terms)[ncol(net_terms)] = condition_name
      rm(condition_pasum)
    }
    rm(condition_idx, condition_name)

    # Select genes for network
    net_genes = NULL
    for (condition_idx in 1:length(network_condition_names)) {
      condition_name = network_condition_names[condition_idx]
      condition_entities = PA_results[[which(names(PA_results) == condition_name)]]$PA$tested_entities[[gene_depth]]
      net_genes = unique(rbind(net_genes, condition_entities))
      rm(condition_name, condition_entities)
    }
    rm(condition_idx)

    colnames(net_genes)[1] = 'UNIPROT'
    net_genes = net_genes[which(!is.na(net_genes$ENTREZID)), ]
    entrez_symbol_map = map_entrez_to_symbols_uniprot(net_genes$ENTREZID, annotation_maps)[, 1:2]
    net_genes = dplyr::left_join(net_genes, entrez_symbol_map, by = "ENTREZID")
    rm(entrez_symbol_map)

    # Mark which conditions each gene appears in
    for (condition_idx in 1:length(network_condition_names)) {
      condition_name = network_condition_names[condition_idx]
      condition_entities = PA_results[[which(names(PA_results) == condition_name)]]$PA$tested_entities[[gene_depth]]$entity
      net_genes$condition_flag = 0
      net_genes$condition_flag[which(net_genes$UNIPROT %in% condition_entities)] = 1
      colnames(net_genes)[ncol(net_genes)] = condition_name
      rm(condition_name, condition_entities)
    }
    rm(condition_idx)

    # Subset term-gene matrix to selected terms and genes
    tg_net = term_gene_matrix[which(rownames(term_gene_matrix) %in% net_genes$ENTREZID), ]
    tg_net = tg_net[, which(colnames(tg_net) %in% net_terms$term.id)]
    tg_net = tg_net[which(Matrix::rowSums(tg_net) > 0), ]

    range(Matrix::colSums(tg_net))
    range(Matrix::rowSums(tg_net))

    tg_net = as.matrix(tg_net)

    net_terms_use = net_terms[which(net_terms$term.id %in% colnames(tg_net)), ]
    net_genes_use = net_genes[which(net_genes$ENTREZID %in% rownames(tg_net)), ]

    print(paste(nrow(net_terms_use), 'terms and', nrow(net_genes_use), 'genes selected for network'))

    # --------------------------------------------------------------------------
    # Network construction
    # --------------------------------------------------------------------------

    # Build edge list
    network_edges = NULL
    for (term_col_idx in which(colSums(tg_net) > 0)) {
      gene_row_indices = which(tg_net[, term_col_idx] > 0)
      term_id = colnames(tg_net)[term_col_idx]
      term_edges = data.frame(from = rownames(tg_net)[gene_row_indices], to = term_id)
      network_edges = rbind(network_edges, term_edges)
      rm(gene_row_indices, term_id, term_edges)
    }
    rm(term_col_idx)

    # Build node list
    network_nodes = data.frame(
      name = c(rownames(tg_net), colnames(tg_net)),
      type = c(rep('protein', nrow(tg_net)), rep('term', ncol(tg_net))),
      label = c(rep('', nrow(tg_net)), colnames(tg_net))
    )

    # Initial color assignment
    network_nodes$color = adjustcolor('black', alpha.f = 0.2)
    network_nodes$label.color = 'black'

    # Identify nodes by condition
    idx_pre_up = which(
      network_nodes$name %in% net_genes$ENTREZID[which(net_genes[[network_condition_names[1]]] == 1)] |
      network_nodes$name %in% net_terms$term.id[which(net_terms[[network_condition_names[1]]] == 1)]
    )
    idx_pre_down = which(
      network_nodes$name %in% net_genes$ENTREZID[which(net_genes[[network_condition_names[2]]] == 1)] |
      network_nodes$name %in% net_terms$term.id[which(net_terms[[network_condition_names[2]]] == 1)]
    )
    idx_post_up = which(
      network_nodes$name %in% net_genes$ENTREZID[which(net_genes[[network_condition_names[3]]] == 1)] |
      network_nodes$name %in% net_terms$term.id[which(net_terms[[network_condition_names[3]]] == 1)]
    )
    idx_post_down = which(
      network_nodes$name %in% net_genes$ENTREZID[which(net_genes[[network_condition_names[4]]] == 1)] |
      network_nodes$name %in% net_terms$term.id[which(net_terms[[network_condition_names[4]]] == 1)]
    )

    # Color by condition (old scheme - will be replaced)
    network_nodes$label.color[idx_pre_up] = 'red'
    network_nodes$label.color[idx_pre_down] = 'blue'
    network_nodes$label.color[intersect(idx_pre_up, idx_pre_down)] = 'purple'

    network_nodes$color[idx_post_up] = adjustcolor('red', alpha.f = 0.5)
    network_nodes$color[idx_post_down] = adjustcolor('blue', alpha.f = 0.5)
    network_nodes$color[intersect(idx_post_up, idx_post_down)] = 'purple'

    # Visual properties
    network_nodes$label.font = 2
    network_nodes$label.font[which(network_nodes$type == 'term')] = 3
    network_nodes$label.family = "Tahoma"
    network_nodes$label.family[which(network_nodes$type == 'term')] = "Arial Narrow, Tahoma"
    network_nodes$size = 5
    network_nodes$size[which(network_nodes$type == 'term')] = 15
    network_nodes$label.cex = 0.3
    network_nodes$label.cex[which(network_nodes$type == 'term')] = 0.5
    network_nodes$shape = 'circle'
    network_nodes$shape[which(network_nodes$type == 'term')] = 'rectangle'

    # Recolor with comprehensive condition scheme
    node_colors = rep('grey60', nrow(network_nodes))
    color_legend = NULL

    color_legend = c(color_legend, gplots::col2hex('palevioletred2'))
    names(color_legend)[length(color_legend)] = 'preDx up'
    node_colors[setdiff(idx_pre_up, c(idx_pre_down, idx_post_up, idx_post_down))] = gplots::col2hex('palevioletred2')

    color_legend = c(color_legend, gplots::col2hex('firebrick2'))
    names(color_legend)[length(color_legend)] = 'preDx up & postDx up'
    node_colors[setdiff(intersect(idx_pre_up, idx_post_up), c(idx_pre_down, idx_post_down))] = gplots::col2hex('firebrick2')

    color_legend = c(color_legend, gplots::col2hex('darkred'))
    names(color_legend)[length(color_legend)] = 'postDx up'
    node_colors[setdiff(idx_post_up, c(idx_pre_up, idx_pre_down, idx_post_down))] = gplots::col2hex('darkred')

    color_legend = c(color_legend, gplots::col2hex('deepskyblue2'))
    names(color_legend)[length(color_legend)] = 'preDx down'
    node_colors[setdiff(idx_pre_down, c(idx_pre_up, idx_post_up, idx_post_down))] = gplots::col2hex('deepskyblue2')

    color_legend = c(color_legend, gplots::col2hex('royalblue'))
    names(color_legend)[length(color_legend)] = 'preDx down & postDx down'
    node_colors[setdiff(intersect(idx_pre_down, idx_post_down), c(idx_pre_up, idx_post_up))] = gplots::col2hex('royalblue')

    color_legend = c(color_legend, gplots::col2hex('navyblue'))
    names(color_legend)[length(color_legend)] = 'postDx down'
    node_colors[setdiff(idx_post_down, c(idx_pre_up, idx_pre_down, idx_post_up))] = gplots::col2hex('navyblue')

    color_legend = c(color_legend, gplots::col2hex('orange'))
    names(color_legend)[length(color_legend)] = 'preDx up & postDx down'
    node_colors[setdiff(intersect(idx_pre_up, idx_post_down), c(idx_pre_down, idx_post_up))] = gplots::col2hex('orange')

    color_legend = c(color_legend, gplots::col2hex('forestgreen'))
    names(color_legend)[length(color_legend)] = 'preDx down & postDx up'
    node_colors[setdiff(intersect(idx_pre_down, idx_post_up), c(idx_pre_up, idx_post_down))] = gplots::col2hex('forestgreen')

    network_nodes$label.color = node_colors
    network_nodes$color = node_colors

    # Assign human-readable labels
    for (node_idx in 1:nrow(network_nodes)) {
      gene_match = which(net_genes$ENTREZID == network_nodes$name[node_idx])
      if (length(gene_match) > 0) {
        network_nodes$label[node_idx] = net_genes$SYMBOL[gene_match]
      }
      term_match = which(net_terms$term.id == network_nodes$name[node_idx])
      if (length(term_match) > 0) {
        network_nodes$label[node_idx] = net_terms$Description[term_match]
      }
      rm(gene_match, term_match)
    }
    rm(node_idx)

    # Build igraph object
    network_graph = igraph::graph_from_data_frame(network_edges, vertices = network_nodes, directed = FALSE)

    print(paste('network:', nrow(network_nodes), 'nodes,', nrow(network_edges), 'edges.'))

    # --------------------------------------------------------------------------
    # Interactive network visualizations with visNetwork
    # --------------------------------------------------------------------------

    # Version 1: No grouping
    igraph_to_visnetwork(network_graph)

    # Version 2: Group proteins with identical term association patterns
    network_graph_v2 = collapse_terms_by_identical_proteins(
      network_graph,
      term_value = 'protein',
      protein_value = 'term',
      representative = 'concat'
    )
    igraph_to_visnetwork(network_graph_v2)

    # Version 3: Two-fold grouping (terms first, then proteins)
    network_graph_v3 = collapse_terms_by_identical_proteins(network_graph)
    network_graph_v3 = collapse_terms_by_identical_proteins(
      network_graph_v3,
      term_value = 'protein',
      protein_value = 'term',
      representative = 'concat'
    )
    igraph_to_visnetwork(network_graph_v3)

    # Plot color legend
    pdf(OUTPUT_FILE_NETWORK_STATIC_LEGEND, width = 7, height = 7)
    plot_kololab_legend(color_legend, title = 'regulation')
    dev.off()

    # --------------------------------------------------------------------------
    # Export static network plot as PDF
    # --------------------------------------------------------------------------
    if (EXPORT_NETWORK_STATIC) {

      names(grDevices::pdfFonts())

      network_graph_export = network_graph
      igraph::V(network_graph_export)$size[which(igraph::V(network_graph_export)$shape == 'circle')] = 4
      igraph::V(network_graph_export)$shape[which(igraph::V(network_graph_export)$shape == 'rectangle')] = 'square'
      igraph::V(network_graph_export)$size[which(igraph::V(network_graph_export)$shape == 'square')] = 6
      igraph::V(network_graph_export)$label.family[which(igraph::V(network_graph_export)$shape == 'circle')] = 'serif'
      igraph::V(network_graph_export)$label.family[which(igraph::V(network_graph_export)$shape == 'square')] = 'sans'

      pdf(OUTPUT_FILE_NETWORK_STATIC, width = 7, height = 7)
      set.seed(1)
      plot(
        network_graph_export,
        vertex.label.dist = -1,
        layout = igraph::layout_with_fr(network_graph_export, niter = 100000)
      )
      dev.off()

      print(paste('Exported:', OUTPUT_FILE_NETWORK_STATIC))

      rm(network_graph_export)
    }

  }  # End if figure_type == 'main'

}  # End for loop (main & supplementary)

}  # END PART 3 (Network generation is now integrated into the for loop above)
