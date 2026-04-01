## SPDX-FileCopyrightText: 2026 Alkahest, Inc. <chronospd_data_access@alkahest.com>
##
## SPDX-License-Identifier: MIT

plot_pair <- function(tbl, symmetric.y = FALSE) {
  required_cols <- c(
    "assay_id_numerator",
    "assay_id_denominator",
    "gene_symbol_numerator",
    "gene_symbol_denominator",
    "ratio",
    "NPX.numerator",
    "NPX.denominator",
    "YEARS_SINCE_ONSET",
    "group"
  )
  missing_cols <- setdiff(required_cols, names(tbl))
  
  if (length(missing_cols) > 0L) {
    stop(
      sprintf("Missing required column(s): %s", paste(missing_cols, collapse = ", ")),
      call. = FALSE
    )
  }

  .d <- tbl |>
    dplyr::filter(
      !is.na(NPX.numerator),
      !is.na(NPX.denominator)
    )
  
  if(nrow(.d) == 0L) {
    return(invisible(NULL))
  }
  
  #####
  # Determine y-axis range(s)
  #####
  #(.d |> dplyr::pull(NPX.numerator) |> range(na.rm = TRUE) |> abs() |> max()) * c(-1, 1)
  
  numerator_plot <- .d |>
    plot_protein(y_var = "NPX.numerator", symmetric.y = symmetric.y) +
    ggplot2::facet_wrap(
      glue::glue("{gene_symbol_numerator}<br /><span style='font-size:7pt;color:gray;'>{assay_id_numerator}</span>") ~ .,
      ncol = 1L,
      scales = "free_y"
    ) +
    ggplot2::theme(
      strip.text.x = ggtext::element_markdown()
    )
  
  denominator_plot <- .d |>
    plot_protein(y_var = "NPX.denominator", symmetric.y = symmetric.y) +
    ggplot2::facet_wrap(
      glue::glue("{gene_symbol_denominator}<br /><span style='font-size:7pt;color:gray;'>{assay_id_denominator}</span>") ~ .,
      #glue::glue("{gene_symbol_denominator}\n{assay_id_denominator}") ~ .,
      ncol = 1L,
      scales = "free_y"
    ) +
    ggplot2::theme(
      strip.text.x = ggtext::element_markdown()
    )
  
  ratio_plots <- .d |>
    plot_protein(y_var = "ratio", y_label = "NPX Ratio", symmetric.y = symmetric.y) +
    ggplot2::facet_wrap(
      glue::glue("log2({gene_symbol_numerator}) - log2({gene_symbol_denominator})<br /><span style='font-size:7pt;color:gray;'>log2({assay_id_numerator}) - log2({assay_id_denominator})</span>") ~ .,
      #glue::glue("log2({gene_symbol_numerator}) - log2({gene_symbol_denominator})\nlog2({assay_id_numerator}) - log2({assay_id_denominator})") ~ .,
      ncol = 1L,
      scales = "free_y"
    ) +
    ggplot2::theme(
      strip.text.x = ggtext::element_markdown()
    )
  
  #####
  # Plot Pair
  #####
  p <- numerator_plot + denominator_plot + ratio_plots +
    patchwork::plot_layout(
      design = "123",
      guides = "collect",
      axis_titles = "collect",
      axes = "collect"
    ) &
    ggplot2::theme(
      legend.position = "bottom",
      legend.title = ggplot2::element_blank()
    )
  
  print(p)
}
