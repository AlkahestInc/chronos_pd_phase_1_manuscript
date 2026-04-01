## SPDX-FileCopyrightText: 2026 Alkahest, Inc. <chronospd_data_access@alkahest.com>
##
## SPDX-License-Identifier: MIT

#' Create volcano plot for differential expression analysis
#'
#' @param dataset_result Processed data frame with columns: gene_symbol, uniprot_id, log2_fc_all, p_val_all, q_val_all, point_group
#' @param paper_pd_list A data frame containing PD hit proteins with uniprot_id column
#' @param x_label Label for x-axis (default: "log2FC between PD and control")
#' @param y_label Label for y-axis (default: "-log10(p.value)")
#' @param x_limits Limits for x-axis (default: c(-0.5, 0.8))
#' @param y_limits Limits for y-axis (default: c(0, 11))
#' @param point_size Size of points (default: 0.7)
#' @param text_size Size of text labels (default: 2.5)
#' @param axis_text_size Size of axis text (default: 15)
#' @param axis_title_size Size of axis title (default: 15)
#' @param max_overlaps Maximum number of overlapping text labels (default: 10)
#' @param color_values Named vector of colors for point groups (default: standard colors)
#'
#' @return A ggplot2 object
create_volcano_plot <- function(
    dataset_result,
    paper_pd_list,
    x_label = "log2FC between PD and control",
    y_label = "-log10(p.value)",
    x_limits = c(-0.5, 0.8),
    y_limits = c(0, 11),
    point_size = 0.7,
    text_size = 2.5,
    axis_text_size = 15,
    axis_title_size = 15,
    max_overlaps = 10,
    color_values = c(
      NS = "grey",
      P_sig = "grey60",
      FDR_sig = "grey20",
      PD_hit = "#FF0066"
    )
) {
  p <- dataset_result |>
    dplyr::mutate(
      point_group = factor(
        point_group,
        levels = c("NS", "P_sig", "FDR_sig", "PD_hit")
      )
    ) |>
    dplyr::arrange(point_group) |>
    ggplot2::ggplot(
      ggplot2::aes(
        x = log2_fc_all,
        y = -log10(p_val_all),
        colour = point_group
      )
    ) +
    ggplot2::geom_point(size = point_size) +
    ggplot2::scale_colour_manual(values = color_values) +
    ggplot2::labs(
      x = x_label,
      y = y_label
    ) +
    ggplot2::scale_y_continuous(limits = y_limits) +
    ggplot2::scale_x_continuous(limits = x_limits) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text = ggplot2::element_text(size = axis_text_size),
      axis.title = ggplot2::element_text(size = axis_title_size, face = "bold"),
      legend.position = "none"
    )

  # Add identifying labels to datapoints
  p <- p +
    ggrepel::geom_text_repel(
      ggplot2::aes(
        label = ifelse(
          uniprot_id %in% paper_pd_list$uniprot_id & q_val_all < 0.05,
          gene_symbol,
          ""
        )
      ),
      col = "red",
      size = text_size,
      max.overlaps = max_overlaps,
      segment.color = NA
    ) +
    ggrepel::geom_text_repel(
      ggplot2::aes(
        label = ifelse(
          (!uniprot_id %in% paper_pd_list$uniprot_id) & q_val_all < 0.01,
          gene_symbol,
          ""
        )
      ),
      col = "black",
      size = text_size,
      max.overlaps = max_overlaps,
      segment.color = NA
    )

  return(p)
}
