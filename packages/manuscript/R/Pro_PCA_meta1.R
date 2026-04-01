## SPDX-FileCopyrightText: 2026 Alkahest, Inc. <chronospd_data_access@alkahest.com>
##
## SPDX-License-Identifier: MIT

Pro_PCA_meta1 <- function(dataset, protein_names, var_names, title_name, cont, nsd, IDs_for_outliers) {
  print(paste("meta examined", var_names))

  Pro_PCA <- prcomp(dataset[, protein_names], center = TRUE, scale = TRUE)
  print(var_names)
  outliers <- apply(Pro_PCA$x[, 1:2], 2, function(x) which(abs(x - median(x)) > (nsd * mad(x)))) ### detection of outliers based on Somalogic's suggestion

  index <- union(unlist(outliers[1]), unlist(outliers[2]))

  if (length(outliers) == 0) {
    print("no outliers")
  } else if (length(outliers) > 0) {
    print(dataset[index, IDs_for_outliers])
  }


  xscore <- as.data.frame(Pro_PCA$x)

  xscore$var <- dataset[, var_names]
  print(paste("dimensions:", dim(xscore)))

  if (cont == "continuous") {
    p <- ggplot2::ggplot(data = xscore, ggplot2::aes(x = PC1, y = PC2, color = var)) +
      ggplot2::geom_point(size = 2, ggplot2::aes(x = PC1, y = PC2)) +
      ggplot2::stat_ellipse(type = "norm", linetype = 2, level = 0.95, ggplot2::aes(x = PC1, y = PC2)) +
      ggplot2::theme_minimal() +
      ggplot2::xlab(paste("Dim1", "(", round(Pro_PCA$sdev[1]^2 / sum(Pro_PCA$sdev^2) * 100, 1), "%", ")", sep = "")) +
      ggplot2::ylab(paste("Dim2", "(", round(Pro_PCA$sdev[2]^2 / sum(Pro_PCA$sdev^2) * 100, 1), "%", ")", sep = "")) +
      ggplot2::labs(col = var_names) +
      ggplot2::ggtitle(paste("2D PCA-plot", "of", title_name)) +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) +
      ggplot2::theme(axis.text = ggplot2::element_text(size = 15), axis.title = ggplot2::element_text(size = 20, face = "bold"))

    print(p)
  } else if (cont != "continuous") {
    p <- ggplot2::ggplot(data = xscore, ggplot2::aes(x = PC1, y = PC2)) +
      ggplot2::geom_point(size = 2, ggplot2::aes(x = PC1, y = PC2, color = as.factor(var))) +
      ggplot2::stat_ellipse(type = "norm", linetype = 2, level = 0.95, ggplot2::aes(x = PC1, y = PC2, color = as.factor(var))) +
      ggplot2::theme_minimal() +
      ggplot2::xlab(paste("Dim1", "(", round(Pro_PCA$sdev[1]^2 / sum(Pro_PCA$sdev^2) * 100, 1), "%", ")", sep = "")) +
      ggplot2::theme(legend.position = "right") +
      ggplot2::ylab(paste("Dim2", "(", round(Pro_PCA$sdev[2]^2 / sum(Pro_PCA$sdev^2) * 100, 1), "%", ")", sep = "")) +
      ggplot2::labs(col = var_names) +
      ggplot2::ggtitle(paste("2D PCA-plot", "of", title_name)) +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) +
      ggplot2::theme(axis.text = ggplot2::element_text(size = 15), axis.title = ggplot2::element_text(size = 20, face = "bold")) +
      ggplot2::scale_color_manual(values = c("purple4", "yellow4"))

    print(p)
  }
}
