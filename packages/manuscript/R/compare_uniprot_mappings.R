## SPDX-FileCopyrightText: 2026 Alkahest, Inc. <chronospd_data_access@alkahest.com>
##
## SPDX-License-Identifier: MIT

# ------------------------------------------------------------------------------
# Function: compare_uniprot_mappings
# Description: Safeguard to check if all PA runs use the same biological entity to ENTREZID mapping
# Parameters:
#   PAs - list of PA run results
# Returns: data frame with pairwise comparison results
# ------------------------------------------------------------------------------
compare_uniprot_mappings <- function(PAs) {
  
  # Pull background tables from PA runs
  extract_background = function(pa_result) {
    background = tryCatch(pa_result$PA$dataset_background, error = function(e) NULL)
    if (is.null(background) || !"entity" %in% names(background) || !"ENTREZID" %in% names(background)) {
      return(NULL)
    }
    background = unique(background[, c("entity", "ENTREZID")])
    background[] = lapply(background, as.character)
    background
  }
  
  background_list = lapply(PAs, extract_background)
  names(background_list) = if (is.null(names(PAs))) {
    paste0("PA", seq_along(PAs))
  } else {
    names(PAs)
  }
  
  # Keep only PA runs that have a background table
  has_background = !vapply(background_list, is.null, logical(1))
  if (!any(has_background)) {
    stop("None of the elements in 'PAs' contains PA$dataset_background.")
  }
  
  skipped = names(background_list)[!has_background]
  background_list = background_list[has_background]
  
  # All possible pairwise comparisons
  pair_indices = utils::combn(seq_along(background_list), 2, simplify = FALSE)
  
  comparison_df = do.call(
    rbind,
    lapply(pair_indices, function(idx) {
      i = idx[1]
      j = idx[2]
      bg1 = background_list[[i]]
      bg2 = background_list[[j]]
      
      shared_entities = intersect(bg1$entity, bg2$entity)
      
      # For shared entities, extract ENTREZID used in PA run
      entrez1 = bg1[match(shared_entities, bg1$entity), "ENTREZID", drop = TRUE]
      entrez2 = bg2[match(shared_entities, bg2$entity), "ENTREZID", drop = TRUE]
      
      data.frame(
        dataset1 = names(background_list)[i],
        dataset2 = names(background_list)[j],
        n_mappings_1 = nrow(bg1),
        n_mappings_2 = nrow(bg2),
        n_shared = length(shared_entities),
        same_size = nrow(bg1) == nrow(bg2),
        same_mapping_in_shared = identical(entrez1, entrez2),
        stringsAsFactors = FALSE
      )
    })
  )
  
  attr(comparison_df, "skipped") = skipped
  comparison_df
}
