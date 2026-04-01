## SPDX-FileCopyrightText: 2026 Alkahest, Inc. <chronospd_data_access@alkahest.com>
##
## SPDX-License-Identifier: MIT

ukb_rap_upload <- function(
  from,
  rap_path,
  project = Sys.getenv("DX_PROJECT_CONTEXT_ID")
) {
  # Ensure the folder already exists
  folder_record_id <- system(
    glue::glue(
      'dx mkdir --parents "{folder}"',
      folder = dirname(rap_path)
    ),
    intern = TRUE
  )

  # Check if the file already exists
  file_list <- tryCatch(
    {
      system(
        glue::glue('dx ls --brief "{rap_path}"'),
        intern = TRUE,
        ignore.stderr = TRUE
      )
    },
    warning = function(w) {
      #message("Warning: File does not exist")
      return(character(0)) # Return an empty character vector
    },
    error = function(e) {
      #message("Error: Unable to resolve file")
      return(character(0)) # Return an empty character vector
    }
  )

  if (length(file_list)) {
    stop(glue::glue("File already exists, won't overwrite it: {rap_path}"))
  }

  result <- system(
    glue::glue(
      'dx upload "{from}" --destination "{project}:{rap_path}"',
      from = from,
      project = project,
      rap_path = rap_path
    ),
    intern = TRUE
  ) |>
    tibble::as_tibble()

  invisible(result)
}
