## SPDX-FileCopyrightText: 2026 Alkahest, Inc. <chronospd_data_access@alkahest.com>
##
## SPDX-License-Identifier: MIT

pick_single_sample_per_participant <- function(.data, strategy = c("middle", "first", "last", "random")) {
  strategy <- match.arg(strategy)

  required_cols <- c(
    "Participant_ID",
    "YEARS_SINCE_ONSET"
  )
  missing_cols <- setdiff(required_cols, names(.data))

  if (length(missing_cols) > 0L) {
    stop(
      sprintf("Missing required column(s): %s", paste(missing_cols, collapse = ", ")),
      call. = FALSE
    )
  }

  if(nrow(.data) <= 0) {
    return(.data)
  }
  
  .data <- .data |>
    dplyr::group_by(
      Participant_ID
    ) |>
    dplyr::arrange(
      YEARS_SINCE_ONSET
    ) |>
    dplyr::mutate(
      selected_sample = dplyr::case_when(
        strategy == "last" ~ dplyr::n(),
        strategy == "first" ~ 1,
        strategy == "random" ~ sample(dplyr::row_number(), 1),
        strategy == "middle" ~ {
          mid_point <- (min(YEARS_SINCE_ONSET) + max(YEARS_SINCE_ONSET)) / 2
          which.min(abs(YEARS_SINCE_ONSET - mid_point))
        },
        TRUE ~ NA_integer_
      )
    ) |>
    dplyr::filter(selected_sample == dplyr::row_number()) |>
    dplyr::select(-selected_sample) |>
    dplyr::ungroup()

  return(.data)
}
