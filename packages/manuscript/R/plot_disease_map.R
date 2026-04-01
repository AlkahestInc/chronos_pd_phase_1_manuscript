## SPDX-FileCopyrightText: 2026 Alkahest, Inc. <chronospd_data_access@alkahest.com>
##
## SPDX-License-Identifier: MIT

plot_disease_map <- function(
  tbl,
  group_theme,
  y = "n_donors",
  y_tick_values = 10L^c(2L:7L)
) {
  group_theme <- group_theme |>
    dplyr::mutate(
      group = factor(group, levels = unique(group)),
      text_x = purrr::map(
        title,
        ~ ggplot2::element_text(colour = .x, face = "bold")
      ),
    )

  tbl <- tbl |>
    dplyr::group_by(group) |>
    dplyr::mutate(
      x = dplyr::row_number()
    ) |>
    dplyr::ungroup()

  p <- tbl |>
    ggplot2::ggplot(
      ggplot2::aes(
        x = x,
        y = .data[[y]],
        label = label,
        color = group
      )
    ) +
    ggplot2::labs(
      title = "Disease Map",
      x = "Indication Spaces",
      y = "Number of Donors",
    ) +
    ggh4x::facet_grid2(
      . ~ group,
      strip = ggh4x::strip_themed(
        text_x = group_theme[["text_x"]]
      ),
      scales = "free_x",
      space = "free_x"
    ) +
    ggplot2::scale_x_continuous(
      expand = ggplot2::expansion(
        mult = c(0, 0),
        # Add fixed padding to both sides of the x-axis
        add = c(1.5, 1.5)
      )
    ) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(
        size = 14,
        face = "bold",
        vjust = 0,
        margin = ggplot2::margin(b = 10)
      ),
      legend.position = "none",
      axis.line = ggplot2::element_line(linewidth = 1.5, color = "black"),
      axis.ticks.y = ggplot2::element_line(linewidth = 1.5, color = "black"),
      axis.ticks.length = grid::unit(10, "pt"),
      strip.background = ggplot2::element_blank(),
      panel.grid.minor.y = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_blank(),
      panel.grid.minor.x = ggplot2::element_blank(),
      panel.spacing = grid::unit(0, "pt"),
      plot.background = ggplot2::element_blank(),
      panel.background = ggplot2::element_rect(fill = "white")
    ) +
    #####
    # Add background colour to each group
    #####
    ggplot2::geom_rect(
      data = group_theme,
      inherit.aes = FALSE,
      ggplot2::aes(
        fill = group
      ),
      linewidth = 0,
      xmin = -Inf,
      xmax = Inf,
      ymin = 0,
      ymax = 10
    ) +
    ggplot2::scale_fill_manual(
      values = setNames(
        group_theme[["background_colour"]],
        group_theme[["group"]]
      )
    ) +
    #####
    # Place horizontal guide lines at the same places as the y-axis tick marks
    #####
    ggplot2::geom_hline(
      yintercept = y_tick_values,
      color = "gray80",
      linewidth = 1.5,
      linetype = "dashed"
    ) +
    #####
    # Customise the y-axis labels
    #####
    ggplot2::scale_y_log10(
      breaks = y_tick_values,
      limits = range(y_tick_values),
      # labels = scales::comma(y_tick_values),
      labels = glue::glue(
        "<img src='{src}' width='{width}'/><br>{text}",
        src = system.file(glue::glue(
          "img/donors_{y_tick_values}.png",
          y_tick_values = format(y_tick_values, scientific = FALSE, trim = TRUE)
        ), package = "manuscript"),
        width = 10 * log10(y_tick_values),
        text = scales::comma(y_tick_values)
      )
    ) +
    ggplot2::theme(
      axis.text.y = ggtext::element_markdown(
        color = "black",
        size = 11,
        hjust = 0.5
      ),
    ) +
    #####
    # Add data points and their labels
    #####
    ggplot2::geom_point(
      size = 7,
    ) +
    ggplot2::scale_colour_manual(
      values = setNames(
        group_theme[["datapoint"]],
        group_theme[["group"]]
      )
    ) +
    ggrepel::geom_text_repel(
      colour = "black",
      size = 3,
      direction = "y", # only nudge text vertically
      box.padding = 0.5, # prevent label-label overlaps
      point.padding = 2, # move labels away from own data point
      max.time = 2,
      max.iter = 10000,
      max.overlaps = Inf, # show all labels
    )

  return(p)
}
