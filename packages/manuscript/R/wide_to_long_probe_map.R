## SPDX-FileCopyrightText: 2026 Alkahest, Inc. <chronospd_data_access@alkahest.com>
##
## SPDX-License-Identifier: MIT

# ------------------------------------------------------------------------------
# Function: wide_to_long_probe_map
# Description: Convert a wide probe annotation table to long format
# Details:
#   Converts a probe-level wide annotation table back to long format by
#   expanding pipe-separated UniProt and gene symbol fields into one row
#   per mapping.
#   Probe-level result columns are preserved and replicated for each
#   reconstructed mapping. Empty gene symbol slots are interpreted as
#   missing values (`NA`). Positional missing symbols are reconstructed
#   when this is unambiguous, for example from leading, trailing, or
#   repeated separators. The function stops with an error if the symbol
#   mapping cannot be resolved uniquely.
# Parameters:
#   x - A data.frame in wide format containing one row per probe.
#       Must include the probe identifier column, the UniProt column, the
#       gene symbol column, and any additional probe-level result columns.
#   probe_col -  Character scalar. Name of the column containing the
#       probe identifier. This column must be unique in `x`.
#   uniprot_col - Character scalar. Name of the column containing
#       pipe-separated UniProt identifiers.
#   symbol_col - Character scalar. Name of the column containing
#       pipe-separated gene symbols. Empty entries are treated as missing
#       values and returned as `NA`.
#   sep - Character scalar. Separator used to split multi-valued
#       annotation fields. Default is `"|"`.
# Returns:
#   A data.frame in long format with one row per probe-to-UniProt
#   mapping. Probe-level columns are preserved, `probe_col` may appear
#   multiple times, `uniprot_col` contains a single UniProt ID per row,
#   and `symbol_col` contains the corresponding gene symbol or `NA`.
# ------------------------------------------------------------------------------
wide_to_long_probe_map <- function(
  x,
  probe_col = "probe_id",
  uniprot_col = "uniprot_id",
  symbol_col = "gene_symbol",
  sep = "|"
) {
  if (!is.data.frame(x)) {
    stop("`x` must be a data.frame.", call. = FALSE)
  }
  
  needed <- c(probe_col, uniprot_col, symbol_col)
  missing_cols <- setdiff(needed, colnames(x))
  if (length(missing_cols) > 0L) {
    stop(
      sprintf(
        "Missing required column(s): %s",
        paste(missing_cols, collapse = ", ")
      ),
      call. = FALSE
    )
  }
  
  if (nrow(x) == 0L) {
    return(x[0, , drop = FALSE])
  }
  
  if (anyDuplicated(x[[probe_col]]) > 0L) {
    dup_ids <- unique(x[[probe_col]][duplicated(x[[probe_col]])])
    stop(
      sprintf(
        "Input is not valid wide format: `%s` is not unique. Duplicated probe(s): %s",
        probe_col,
        paste(utils::head(as.character(dup_ids), 10L), collapse = ", ")
      ),
      call. = FALSE
    )
  }
  
  x <- as.data.frame(x, stringsAsFactors = FALSE)
  
  # fontos: ezek karakterek legyenek, különben faktor-szinteknél assignment gond lehet
  x[[uniprot_col]] <- as.character(x[[uniprot_col]])
  x[[symbol_col]] <- as.character(x[[symbol_col]])
  
  fmt1 <- function(z) {
    if (length(z) == 0L || is.na(z)) "NA" else paste0("'", z, "'")
  }
  
  trim_and_na <- function(v) {
    v <- trimws(v)
    v[v == ""] <- NA_character_
    v
  }
  
  # split, de az üres elemeket is őrizze meg:
  # "|A|B"   -> "", "A", "B"
  # "A|B|"   -> "A", "B", ""
  # "A||B"   -> "A", "", "B"
  # "||"     -> "", "", ""
  split_preserve_empty <- function(s, sep = "|") {
    if (length(s) != 1L) {
      stop("Internal error: `split_preserve_empty()` expects scalar input.", call. = FALSE)
    }
    
    s <- as.character(s)
    
    if (is.na(s)) {
      return(NA_character_)
    }
    
    if (!nzchar(s)) {
      return("")
    }
    
    pos <- gregexpr(sep, s, fixed = TRUE)[[1]]
    
    if (length(pos) == 1L && pos[1] == -1L) {
      return(s)
    }
    
    sep_n <- nchar(sep, type = "chars")
    starts <- c(1L, pos + sep_n)
    ends <- c(pos - 1L, nchar(s, type = "chars"))
    
    substring(s, starts, ends)
  }
  
  out_list <- lapply(seq_len(nrow(x)), function(i) {
    probe_val <- x[[probe_col]][i]
    
    raw_u <- x[[uniprot_col]][i]
    if (is.na(raw_u) || trimws(raw_u) == "") {
      stop(
        sprintf(
          "Probe %s: `%s` is missing or empty.",
          fmt1(as.character(probe_val)),
          uniprot_col
        ),
        call. = FALSE
      )
    }
    
    u <- trim_and_na(split_preserve_empty(raw_u, sep = sep))
    
    if (length(u) == 0L || anyNA(u)) {
      stop(
        sprintf(
          "Probe %s: invalid `%s` value %s. UniProt slots must be explicit and non-empty after splitting.",
          fmt1(as.character(probe_val)),
          uniprot_col,
          fmt1(raw_u)
        ),
        call. = FALSE
      )
    }
    
    n_u <- length(u)
    
    raw_s <- x[[symbol_col]][i]
    
    # teljesen üres / NA symbol mező => minden symbol hiányzik
    if (is.na(raw_s) || trimws(raw_s) == "") {
      s <- rep(NA_character_, n_u)
    } else {
      s <- trim_and_na(split_preserve_empty(raw_s, sep = sep))
      n_s <- length(s)
      
      if (n_s != n_u) {
        if (n_s < n_u) {
          stop(
            sprintf(
              paste0(
                "Probe %s: parsed %d UniProt ID(s) from %s, but only %d symbol slot(s) from %s.\n",
                "This is ambiguous: at least one symbol is missing, but its position cannot be reconstructed.\n",
                "Use explicit empty slots with '%s', e.g. leading/trailing/double separators."
              ),
              fmt1(as.character(probe_val)),
              n_u,
              fmt1(raw_u),
              n_s,
              fmt1(raw_s),
              sep
            ),
            call. = FALSE
          )
        } else {
          stop(
            sprintf(
              paste0(
                "Probe %s: parsed %d UniProt ID(s) from %s, but %d symbol slot(s) from %s.\n",
                "This is inconsistent: there are more symbol slots than UniProt IDs."
              ),
              fmt1(as.character(probe_val)),
              n_u,
              fmt1(raw_u),
              n_s,
              fmt1(raw_s)
            ),
            call. = FALSE
          )
        }
      }
    }
    
    block <- x[rep(i, n_u), , drop = FALSE]
    block[[uniprot_col]] <- u
    block[[symbol_col]] <- s
    
    # extra sanity check
    if (anyNA(block[[uniprot_col]]) || any(trimws(block[[uniprot_col]]) == "")) {
      stop(
        sprintf(
          "Probe %s: internal sanity check failed; empty UniProt generated in output.",
          fmt1(as.character(probe_val))
        ),
        call. = FALSE
      )
    }
    
    block
  })
  
  out <- do.call(rbind, out_list)
  rownames(out) <- NULL
  out
}
