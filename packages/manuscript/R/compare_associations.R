## SPDX-FileCopyrightText: 2026 Alkahest, Inc. <chronospd_data_access@alkahest.com>
##
## SPDX-License-Identifier: MIT

#' Compare Models with Different Associations
#'
#' @param data joint_data object
#' @param spec Base JointModelSpec
#' @param proteins Proteins to test
#' @param associations List of association types to compare
#' @export
compare_associations <- function(data, spec, proteins,
                                 associations = list("value", "slope", c("value", "slope"))) {
  cli::cli_h1("Comparing Association Structures")

  # Fit models for each association type
  comparison_results <- purrr::map(associations, function(assoc) {
    cli::cli_alert_info("Fitting models with {paste(assoc, collapse = ' + ')}")

    spec_copy <- spec$clone()
    spec_copy$set_association(assoc)

    fit_joint_models(data, spec_copy, proteins, progress = FALSE)
  })

  # Extract comparison metrics
  comparison_df <- purrr::map_dfr(seq_along(associations), function(i) {
    results <- comparison_results[[i]]
    assoc_name <- paste(associations[[i]], collapse = "+")

    results$summaries |>
      dplyr::filter(term %in% c("value", "slope")) |>
      dplyr::select(protein, term, dic, waic, lpml) |>
      dplyr::mutate(association = assoc_name)
  })

  # Create comparison visualization
  p_dic <- ggplot2::ggplot(
    comparison_df,
    ggplot2::aes(x = association, y = dic, group = protein)
  ) +
    ggplot2::geom_line(alpha = 0.3) +
    ggplot2::geom_point() +
    ggplot2::facet_wrap(~protein, scales = "free_y") +
    ggplot2::labs(title = "DIC Comparison", y = "DIC (lower is better)") +
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

  list(
    results = comparison_results,
    metrics = comparison_df,
    plot = p_dic,
    best_models = comparison_df |>
      dplyr::group_by(protein) |>
      dplyr::slice_min(dic, n = 1) |>
      dplyr::select(protein, best_association = association, dic)
  )
}
