## SPDX-FileCopyrightText: 2026 Alkahest, Inc. <chronospd_data_access@alkahest.com>
##
## SPDX-License-Identifier: MIT

ukb_rap_get_cohort_data <- local({
  .ukb_rap_get_cohort_data <- function(
    fields,
    cohort_rap_path
  ) {
    fields_file <- tempfile(fileext = ".fields")
    output_file <- tempfile(fileext = ".csv")

    #####
    # Create a file of field names required for constructing the SQL query
    #####
    fields |>
      dplyr::select(field) |>
      readr::write_tsv(
        fields_file,
        col_names = FALSE
      )

    system(
      glue::glue(
        'dx extract_dataset "{cohort_rap_path}" --fields-file "{fields_file}" --output "{output_file}"'
      ),
      intern = FALSE
    )

    print(glue::glue("dx extract_dataset output: {output_file}"))

    new_name <- fields$label[fields$field == "participant.eid"]

    long_data <- readr::read_csv(
      output_file,
      show_col_types = FALSE
    ) |>
      (\(x) {
        other_cols <- setdiff(names(x), "participant.eid")
        if (length(other_cols) > 0) {
          tidyr::pivot_longer(
            x,
            cols = tidyselect::all_of(other_cols),
            names_to = "variable",
            values_to = "coding"
          ) |>
            tidyr::drop_na(coding) |>
            tidyr::separate(
              col = variable,
              into = c("entity_name", "field_name"),
              sep = "\\.", # dot needs to be escaped
              remove = FALSE # keep original 'variable' column if you want
            ) |>
            tidyr::separate_wider_regex(
              field_name,
              c(
                "p",
                field_id = "\\d+",
                "_i",
                instance = "\\d+",
                "_a",
                array = "\\d+"
              ),
              too_few = "align_start",
              cols_remove = FALSE
            ) |>
            dplyr::left_join(
              fields |> dplyr::select(field, label),
              by = c("variable" = "field")
            ) |>
            dplyr::mutate(
              dplyr::across(c(field_id, instance, array), as.integer),
              dplyr::across(c(variable, entity_name, field_name), as.factor),
              field_id_label = as.factor(stringr::str_extract(
                label,
                "^.+?(?= \\|)"
              ))
            )
        } else {
          x
        }
      })() |>
      dplyr::mutate(
        participant.eid = as.integer(participant.eid)
      ) |>
      dplyr::rename(
        !!new_name := participant.eid,
      )

    return(long_data)
  }

  cached_fun <- memoise::memoise(.ukb_rap_get_cohort_data)

  function(field, cohort_rap_path, use_cache = TRUE) {
    if (!use_cache) {
      # Recreate memoised function to invalidate the cache
      cached_fun <<- memoise::memoise(.ukb_rap_get_cohort_data)
    }

    cached_fun(field, cohort_rap_path)
  }
})
