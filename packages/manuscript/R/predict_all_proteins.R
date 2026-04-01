## SPDX-FileCopyrightText: 2026 Alkahest, Inc. <chronospd_data_access@alkahest.com>
##
## SPDX-License-Identifier: MIT

predict_all_proteins <- function(data, protein_cols, years_pred) {
  # Filter to only proteins that exist in the data
  available_proteins <- intersect(protein_cols, colnames(data))

  if (length(available_proteins) == 0) {
    logger::log_info("No proteins found in data")
    return(NULL)
  }

  missing <- setdiff(protein_cols, available_proteins)
  if (length(missing) > 0) {
    logger::log_info("{length(missing)} proteins not found in this dataset")
    logger::log_trace("missing proteins: {paste(missing, collapse = ', ')}")
  }

  # Process all proteins at once
  results <- available_proteins |>
    purrr::map_dfr(
      ~ {
        data |>
          dplyr::rename(protein_value = tidyselect::all_of(.x)) |>
          dplyr::mutate(scaled_protein = as.vector(scale(protein_value))) |>
          dplyr::select(years_since_pd_onset, group, scaled_protein) |>
          tidyr::drop_na() |>
          dplyr::group_by(group) |>
          tidyr::nest() |>
          dplyr::mutate(
            prediction = purrr::map(
              data,
              ~ {
                .x |>
                  dplyr::arrange(years_since_pd_onset) |>
                  {\(d) {
                    loess(scaled_protein ~ years_since_pd_onset, data = d) |>
                      predict(newdata = years_pred)
                  }}()
              }
            )
          ) |>
          dplyr::mutate(probe_id = .x)
      },
      .id = NULL,
      .progress = TRUE
    )

  plx.tot <- results |>
    tidyr::pivot_wider(
      id_cols = probe_id,
      names_from = group,
      values_from = prediction
    ) |>
    dplyr::mutate(
      difference = purrr::map2(
        pd, control,
        ~ .x - .y
      ),
      years_pred = list(years_pred)
    ) |>
    dplyr::select(probe_id, years_pred, difference)

  return(plx.tot)
}
