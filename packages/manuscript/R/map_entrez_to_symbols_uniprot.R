## SPDX-FileCopyrightText: 2026 Alkahest, Inc. <chronospd_data_access@alkahest.com>
##
## SPDX-License-Identifier: MIT

# ------------------------------------------------------------------------------
# Function: map_entrez_to_symbols_uniprot
# Description: Map ENTREZID to SYMBOL and UNIPROT (PA native is ENTREZ)
# Parameters:
#   entrez_ids - character or numeric vector of ENTREZ IDs
#   tg_maps - output of build_TG_for_PA()
#   collapse - delimiter to collapse multi-maps (default: "|")
# Returns: data frame with ENTREZID, SYMBOL, UNIPROT columns
# ------------------------------------------------------------------------------
map_entrez_to_symbols_uniprot <- function(entrez_ids, tg_maps, collapse = "|") {
  
  ids = as.character(entrez_ids)
  
  format_values = function(value_vector) {
    if (is.null(value_vector) || !length(value_vector)) return(NA_character_)
    value_vector = unique(value_vector[!is.na(value_vector) & nzchar(value_vector)])
    if (!length(value_vector)) return(NA_character_)
    paste(value_vector, collapse = collapse)
  }
  
  symbols = vapply(tg_maps$SYMBOL_map[ids], format_values, character(1L))
  uniprots = vapply(tg_maps$UNIPROT_map[ids], format_values, character(1L))
  
  data.frame(
    ENTREZID = ids,
    SYMBOL = symbols,
    UNIPROT = uniprots,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}
