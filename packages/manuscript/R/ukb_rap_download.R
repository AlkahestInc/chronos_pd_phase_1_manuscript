## SPDX-FileCopyrightText: 2026 Alkahest, Inc. <chronospd_data_access@alkahest.com>
##
## SPDX-License-Identifier: MIT

#' Download data from UKB-RAP
#'
#' This function downloads data from UKB-RAP using the `dx` command-line tool.
#'
#' @param rap_path Character. The remote path for the file within the UKB-RAP project.
#' @param local_path Character. The local destination path for the downloaded file.
#' @param project Character. The UKB-RAP project ID. Defaults to the value of the `DX_PROJECT_CONTEXT_ID` environment variable.
#' @param overwrite Logical. Whether to overwrite the destination file if it already exists. Defaults to `FALSE`.
#'
#' @return An invisible tibble containing the command output.
#' @examples
#' \dontrun{
#' ukb_rap_download("remote/path/to/file", "local/path/to/file")
#' }
#' @export
ukb_rap_download <- function(
  rap_path,
  local_path,
  project = Sys.getenv("DX_PROJECT_CONTEXT_ID"),
  overwrite = FALSE
) {
  # Input validation
  if (!is.character(rap_path) || length(rap_path) != 1) {
    stop("'rap_path' must be a single character string.")
  }
  if (!is.character(local_path) || length(local_path) != 1) {
    stop("'local_path' must be a single character string.")
  }
  if (!is.character(project) || length(project) != 1 || project == "") {
    stop("'project' must be a non-empty character string.")
  }
  if (!is.logical(overwrite) || length(overwrite) != 1) {
    stop("'overwrite' must be a single logical value (TRUE or FALSE).")
  }

  system(
    glue::glue(
      'mkdir --parents "{dir}"',
      dir = dirname(local_path)
    ),
    intern = TRUE
  )

  result <- system(
    glue::glue(
      'dx download "{project}:{rap_path}" --output "{local_path}"{if (overwrite) " --overwrite" else ""}'
    ),
    intern = TRUE
  ) |>
    tibble::as_tibble()

  invisible(result)
}
