## SPDX-FileCopyrightText: 2026 Alkahest, Inc. <chronospd_data_access@alkahest.com>
##
## SPDX-License-Identifier: MIT

plot_protein <- function(tbl, x_var = "YEARS_SINCE_ONSET", y_var = "NPX", group_var = "group", x_label = "Years Since Onset", y_label = "NPX", symmetric.y = FALSE) {
  required_cols <- c(
    x_var,
    y_var,
    group_var
  )
  missing_cols <- setdiff(required_cols, names(tbl))
  
  if (length(missing_cols) > 0L) {
    stop(
      sprintf("Missing required column(s): %s", paste(missing_cols, collapse = ", ")),
      call. = FALSE
    )
  }

  p <- tbl |>
    ggplot2::ggplot(
      ggplot2::aes(
        x = !!rlang::sym(x_var),
        y = !!rlang::sym(y_var), 
        colour = !!rlang::sym(group_var),
      )
    ) +
    ggplot2::geom_point(
      shape = "circle",
      stroke = 0.0,
      alpha = 0.3,
      size = 1.0
    ) +
    ggplot2::geom_smooth(
      method = "loess",
      formula = y ~ x,
      method.args = list(
        span = 0.75
      ),
      se = TRUE,
      linewidth = 0.5
    ) +
    ggplot2::scale_color_manual(
      values = c(
        "case.eid" = "darkorange",
        "control.eid"="dodgerblue"
      ),
      labels = c(
        "case.eid" = "Case",
        "control.eid" = "Control"
      )
    ) +
    ggplot2::labs(
      x = x_label,
      y = y_label
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      legend.position = "bottom"
    ) +
    ggplot2::coord_cartesian(
      xlim = c(NA, 5L)
    )
  
  if(symmetric.y) {
    p <- p +
      ggplot2::expand_limits(
        y = (tbl |> dplyr::pull(!!rlang::sym(y_var)) |> range(na.rm = TRUE) |> abs() |> max()) * c(-1, 1)
      )
  }

  return(p)
}
