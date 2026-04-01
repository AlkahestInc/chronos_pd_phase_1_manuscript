## SPDX-FileCopyrightText: 2026 Alkahest, Inc. <chronospd_data_access@alkahest.com>
##
## SPDX-License-Identifier: MIT

#' Print method for joint_data
#' @export
print.joint_data <- function(x, ...) {
  cli::cli_h1("Joint Model Data")
  cli::cli_li("{.val {x$n_subjects}} subjects")
  cli::cli_li("{.val {x$n_observations}} longitudinal observations")
  cli::cli_li("ID variable: {.field {x$id_var}}")
  cli::cli_li("Time variable: {.field {x$time_var}}")

  if (!is.null(x$formula)) {
    cli::cli_h2("Formulas")
    cli::cli_li("Longitudinal: {.code {deparse(x$formula$longitudinal)}}")
    cli::cli_li("Survival: {.code {deparse(x$formula$survival)}}")
  }
  invisible(x)
}
