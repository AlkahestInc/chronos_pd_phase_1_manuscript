## SPDX-FileCopyrightText: 2026 Alkahest, Inc. <chronospd_data_access@alkahest.com>
##
## SPDX-License-Identifier: MIT

# ------------------------------------------------------------------------------
# Function: prepare_probe_annotations_sparse
# Description: Clean probe annotations to single target per probe using scoring rules
# Parameters:
#   df - data frame with probe annotations
#   probe.id - column name for probe IDs
#   annot.col - column name for annotation IDs (e.g., UNIPROT)
#   maps - output from load_annotation_mappings_sparse() (if NULL, loads fresh)
# Returns: list with cleaned data frame and maps
# ------------------------------------------------------------------------------
prepare_probe_annotations_sparse <- function(df,
                                             probe.id,
                                             annot.col,
                                             maps = NULL) {
  # Store original column order
  original_cols <- names(df)
  
  # For safety
  df <- dplyr::ungroup(df)
  
  if (is.null(maps)) maps <- load_annotation_mappings_sparse()
  score_map <- maps$score_map
  
  # Pull out temporarily under ".probe" / ".annot"
  temp_df <- df |>
    dplyr::mutate(
      .probe = .data[[probe.id]],
      .annot = .data[[annot.col]]
    ) |>
    dplyr::ungroup()
  
  # Reliable single-target probes
  reliable_annots <- temp_df |>
    dplyr::group_by(.probe) |>
    dplyr::filter(dplyr::n() == 1) |>
    dplyr::pull(.annot) |>
    unique()
  
  # Split singles vs multis
  probe_counts <- temp_df |> dplyr::count(.probe)
  single_probes <- temp_df |>
    dplyr::filter(.probe %in% probe_counts$.probe[probe_counts$n == 1])
  multi_probes <- temp_df |>
    dplyr::filter(.probe %in% probe_counts$.probe[probe_counts$n > 1])
  
  # Rule 1: prefer annots covered by single-target probes to avoid false positive bubble
  step1_choices <- multi_probes |>
    dplyr::group_by(.probe) |>
    dplyr::summarise(
      chosen = {
        candidates <- .annot
        good_candidates <- intersect(candidates, reliable_annots)
        if (length(good_candidates) == 1) good_candidates else NA_character_
      },
      .groups = "drop"
    )
  
  # Rule 2: prefer non-NA scores, then larger score (more associated terms), then alphabetical
  to_resolve <- step1_choices |>
    dplyr::filter(is.na(chosen)) |>
    dplyr::pull(.probe)
  
  if (length(to_resolve)) {
    step2_choices <- multi_probes |>
      dplyr::filter(.probe %in% to_resolve) |>
      dplyr::distinct(.probe, .annot) |>
      dplyr::mutate(score = score_map[.annot]) |>
      dplyr::group_by(.probe) |>
      dplyr::arrange(dplyr::desc(!is.na(score)), dplyr::desc(score), .annot, .by_group = TRUE) |>
      dplyr::slice_head(n = 1) |>
      dplyr::ungroup() |>
      dplyr::select(.probe, chosen = .annot)
  } else {
    step2_choices <- tibble::tibble(.probe = character(), chosen = character())
  }
  
  # Combine all choices
  final_map <- dplyr::bind_rows(
    single_probes |> dplyr::transmute(.probe, chosen = .annot),
    step1_choices |> dplyr::filter(!is.na(chosen)),
    step2_choices
  )
  
  # Join the chosen annotation back into df, keep original order
  cleaned_df <- df |>
    dplyr::distinct(dplyr::across(dplyr::all_of(probe.id)), .keep_all = TRUE) |>
    dplyr::select(-dplyr::all_of(annot.col)) |>
    dplyr::inner_join(
      final_map |>
        dplyr::rename(
          !!probe.id := .probe,
          !!annot.col := chosen
        ),
      by = probe.id
    ) |>
    dplyr::select(dplyr::all_of(original_cols))
  
  # Output
  list(
    data = cleaned_df,
    maps = maps
  )
}
