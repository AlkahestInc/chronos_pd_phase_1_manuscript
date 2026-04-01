## SPDX-FileCopyrightText: 2026 Alkahest, Inc. <chronospd_data_access@alkahest.com>
##
## SPDX-License-Identifier: MIT

#' Create Joint Model Data Object
#'
#' Modern constructor for joint model data with validation
#'
#' @param longitudinal Tibble with longitudinal measurements
#' @param survival Tibble with survival outcomes
#' @param formula List with longitudinal and survival formulas
#' @param id_var Symbol or string for subject identifier
#' @param time_var Symbol or string for time variable
#'
#' @return joint_data object with S3 class
#' @export
joint_data <- function(longitudinal, survival, formula = NULL,
                       id_var = "id", time_var = "time") {
  # Modern argument handling
  id_var <- rlang::ensym(id_var)
  time_var <- rlang::ensym(time_var)

  # Validate inputs using tidyverse style
  longitudinal <- longitudinal |>
    dplyr::as_tibble() |>
    dplyr::arrange({{ id_var }}, {{ time_var }})

  survival <- survival |>
    dplyr::as_tibble() |>
    dplyr::distinct({{ id_var }}, .keep_all = TRUE)

  # Check consistency
  long_ids <- dplyr::pull(longitudinal, {{ id_var }}) |> unique()
  surv_ids <- dplyr::pull(survival, {{ id_var }})

  if (!setequal(long_ids, surv_ids)) {
    cli::cli_abort(c(
      "Inconsistent subject IDs between longitudinal and survival data",
      "i" = "Longitudinal has {length(setdiff(long_ids, surv_ids))} unique IDs",
      "i" = "Survival has {length(setdiff(surv_ids, long_ids))} unique IDs"
    ))
  }

  structure(
    list(
      longitudinal = longitudinal,
      survival = survival,
      formula = formula,
      id_var = rlang::as_string(id_var),
      time_var = rlang::as_string(time_var),
      n_subjects = length(surv_ids),
      n_observations = nrow(longitudinal)
    ),
    class = c("joint_data", "list")
  )
}
