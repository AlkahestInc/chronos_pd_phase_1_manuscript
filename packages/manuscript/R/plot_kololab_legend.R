## SPDX-FileCopyrightText: 2026 Alkahest, Inc. <chronospd_data_access@alkahest.com>
##
## SPDX-License-Identifier: MIT

# ------------------------------------------------------------------------------
# Function: plot_kololab_legend
# Description: Plot color legend for graph visualization
# Parameters:
#   color_labels - named color vector (names = labels)
#   title - legend title (default: "my legend")
# Returns: invisible(NULL) - draws legend plot
# ------------------------------------------------------------------------------
plot_kololab_legend <- function(color_labels, title = "my legend") {
  
  if (is.null(names(color_labels))) {
    stop("color_labels must be a *named* color vector (names = labels).")
  }
  
  # Keep first occurrence of each color
  unique_indices = !duplicated(color_labels)
  colors = color_labels[unique_indices]
  labels = names(color_labels)[unique_indices]
  
  original_par = par(no.readonly = TRUE)
  on.exit(par(original_par))
  
  par(mar = c(0, 0, if (nzchar(title)) 2 else 0, 0))
  plot.new()
  
  legend(
    "center",
    legend = labels,
    pch = 21,
    pt.bg = colors,
    pt.cex = 2,
    col = "black",
    cex = 1.1,
    bty = "n",
    title = title
  )
}
