## SPDX-FileCopyrightText: 2026 Alkahest, Inc. <chronospd_data_access@alkahest.com>
##
## SPDX-License-Identifier: MIT

LME_results_P1_plot <- function(dataset, coef_names, variable_names, xlim_min, xlim_max, ylim_min, ylim_max, plot, topN, pvalue) {
  LME_results <- dataset |>
    dplyr::filter(coef.names == {{ coef_names }}) |>
    dplyr::select(
      tidyselect::all_of(variable_names)
    )

  LME_results$adjP <- p.adjust(LME_results$p.value, method = "BH")

  LME_results$symbol <- stringr::str_remove(LME_results$protein, "\\..*$") 
  # LME_results$UniProt=UniProt_name
  print(head(LME_results))

  NP005 <- LME_results |>
    dplyr::filter(p.value < 0.05) |>
    dplyr::count()

  NP001 <- LME_results |>
    dplyr::filter(p.value < 0.01) |>
    dplyr::count()

  NP0005 <- LME_results |>
    dplyr::filter(p.value < 0.0005) |>
    dplyr::count()

  NadjP005 <- LME_results |>
    dplyr::filter(adjP < 0.05) |>
    dplyr::count()

  NadjP015 <- LME_results |>
    dplyr::filter(adjP < 0.15) |>
    dplyr::count()

  print(paste(
    "P<0.05:", NP005, "(", round((NP005 / nrow(LME_results)) * 100, 1), "%", ")",
    "P<0.01:", NP001, "(", round((NP001 / nrow(LME_results)) * 100, 1), "%", ")",
    "P<0.0005:", NP0005, "(", round((NP0005 / nrow(LME_results)) * 100, 1), "%", ")",
    "adjP<0.05:", NadjP005, "adjP<0.15:", NadjP015
  ))

  protein_list <- LME_results |>
    dplyr::filter(p.value < 0.05) |>
    dplyr::select(protein, Value, p.value, adjP, symbol) |>
    dplyr::arrange(p.value)

  print(paste("xlim_min:", round(min(LME_results$Value), 3), "xlim_max:", round(max(LME_results$Value), 3), "ylim_max:", round(max(-log10(LME_results$p.value)), 3)))
  print(protein_list |> head(topN))

  if (plot) {
    p <- ggplot2::ggplot(LME_results, ggplot2::aes(x = Value, y = -log10(p.value))) +
      ggplot2::geom_point() +
      ggplot2::xlim(xlim_min, xlim_max) +
      ggplot2::ylim(ylim_min, ylim_max) +
      ggplot2::xlab(paste("coef of ", coef_names)) +
      ggplot2::theme_minimal() +
      ggplot2::geom_point(data = LME_results |> dplyr::filter(p.value < pvalue), col = "red") +
      ggplot2::theme(axis.text = ggplot2::element_text(size = 20), axis.title = ggplot2::element_text(size = 20, face = "bold")) #+
    # annotate(geom="text", x=xlim_min+0.4*(xlim_max-xlim_min), y=ylim_max-2, size=8,
    # label=paste("P<0.05:",NP005,"(",round((NP005/nrow(LME_results))*100,1),"%",")",sep=""),color="red")

    return(p)
  } else {
    return(protein_list)
  }
}
