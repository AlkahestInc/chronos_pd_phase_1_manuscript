## SPDX-FileCopyrightText: 2026 Alkahest, Inc. <chronospd_data_access@alkahest.com>
##
## SPDX-License-Identifier: MIT

# ------------------------------------------------------------------------------
# Function: igraph_to_visnetwork
# Description: Port igraph objects to visNetwork for interactive visualization
# Parameters:
#   graph - igraph object
#   map_shapes - map igraph shapes to visNetwork (circle->dot, rectangle->square) (default: TRUE)
#   label_position - label position: "auto", "below", "inside" (default: "auto")
#   square_offset_px - vertical offset for square labels (default: -10)
#   dot_offset_px - vertical offset for dot labels (default: -5)
# Returns: visNetwork object
# ------------------------------------------------------------------------------
igraph_to_visnetwork <- function(
  graph,
  map_shapes = TRUE,
  label_position = c("auto", "below", "inside"),
  square_offset_px = -10,
  dot_offset_px = -5
) {
  
  stopifnot(inherits(graph, "igraph"))
  if (!requireNamespace("visNetwork", quietly = TRUE)) {
    stop("need to install 'visNetwork' package!")
  }
  
  label_position = match.arg(label_position)
  
  # ============================================================================
  # Build nodes data frame
  # ============================================================================
  
  nodes = igraph::as_data_frame(graph, what = "vertices")
  node_ids = igraph::V(graph)$name
  nodes$id = node_ids
  vertex_attr_names = igraph::vertex_attr_names(graph)
  
  # Node label fallback to id
  if (!"label" %in% names(nodes)) nodes$label = nodes$id
  
  # Map igraph label.font to HTML tags (1=plain, 2=bold, 3=italic, 4=bold+italic)
  label_font = if ("label.font" %in% vertex_attr_names) igraph::V(graph)$label.font else NULL
  if (!is.null(label_font)) {
    for (i in seq_len(nrow(nodes))) {
      font_index = as.integer(label_font[i])
      label_text = nodes$label[i]
      if (is.na(font_index) || is.na(label_text)) next
      
      if (font_index == 2L) {
        label_text = paste0("<b>", label_text, "</b>")
      } else if (font_index == 3L) {
        label_text = paste0("<i>", label_text, "</i>")
      } else if (font_index == 4L) {
        label_text = paste0("<b><i>", label_text, "</i></b>")
      }
      nodes$label[i] = label_text
    }
  }
  
  # Get vertex size
  if ("size" %in% vertex_attr_names) {
    nodes$size = igraph::V(graph)$size
  } else {
    nodes$size = 20
  }
  
  # Shape mapping
  if (!"shape" %in% names(nodes)) nodes$shape = "circle"
  if (map_shapes) {
    shape_values = tolower(as.character(nodes$shape))
    shape_values[shape_values == "circle"] = "dot"
    shape_values[shape_values == "rectangle"] = "square"
    nodes$shape = shape_values
  }
  
  # Translate vertex colors
  fill_color = if ("color" %in% vertex_attr_names) igraph::V(graph)$color else NULL
  fill_color2 = if ("fillcolor" %in% vertex_attr_names) igraph::V(graph)$fillcolor else NULL
  border_color = if ("frame.color" %in% vertex_attr_names) igraph::V(graph)$frame.color else NULL
  
  if (!is.null(fill_color) || !is.null(fill_color2)) {
    background = if (!is.null(fill_color)) fill_color else fill_color2
    if (!is.null(border_color) && length(border_color) == nrow(nodes)) {
      nodes$color = lapply(seq_len(nrow(nodes)), function(i) {
        list(background = background[i], border = border_color[i])
      })
    } else {
      nodes$color = background
    }
  }
  
  # Translate label fonts
  nodes$font = vector("list", nrow(nodes))
  label_color = if ("label.color" %in% vertex_attr_names) igraph::V(graph)$label.color else NULL
  label_cex = if ("label.cex" %in% vertex_attr_names) igraph::V(graph)$label.cex else NULL
  label_family = if ("label.family" %in% vertex_attr_names) igraph::V(graph)$label.family else NULL
  
  # Nudge label position
  is_square = tolower(nodes$shape) == "square"
  is_dot = tolower(nodes$shape) == "dot"
  vertical_adjust = rep(NA_integer_, nrow(nodes))
  
  if (label_position == "auto") {
    vertical_adjust[is_dot] = as.integer(dot_offset_px)
    vertical_adjust[is_square] = as.integer(square_offset_px)
  } else if (label_position == "below") {
    base_offset = pmax(1, as.integer(round(nodes$size * 0.4)))
    vertical_adjust = base_offset
    if (!is.null(dot_offset_px)) vertical_adjust[is_dot] = base_offset[is_dot] + as.integer(dot_offset_px)
    if (!is.na(square_offset_px)) vertical_adjust[is_square] = base_offset[is_square] + as.integer(square_offset_px)
  } else {
    vertical_adjust[] = 0
  }
  
  # Set font properties
  for (i in seq_len(nrow(nodes))) {
    font_props = list()
    if (!is.null(label_color) && !is.na(label_color[i])) font_props$color = label_color[i]
    if (!is.null(label_cex) && !is.na(label_cex[i])) font_props$size = round(label_cex[i] * 48)
    if (!is.null(label_family) && !is.na(label_family[i])) font_props$face = label_family[i]
    if (!is.na(vertical_adjust[i])) font_props$vadjust = vertical_adjust[i]
    if (!is.null(label_font)) font_props$multi = "html"
    nodes$font[[i]] = if (length(font_props)) font_props else NULL
  }
  
  # ============================================================================
  # Build edges data frame
  # ============================================================================
  
  edges = igraph::as_data_frame(graph, what = "edges")
  edge_attr_names = igraph::edge_attr_names(graph)
  
  if (nrow(edges) > 0) {
    if ("color" %in% edge_attr_names) {
      color_values = igraph::E(graph)$color
      if (length(color_values) == nrow(edges)) edges$color = color_values
    }
    if ("width" %in% edge_attr_names) {
      width_values = igraph::E(graph)$width
      if (length(width_values) == nrow(edges)) edges$width = width_values
    }
    if ("label" %in% edge_attr_names) {
      label_values = igraph::E(graph)$label
      if (length(label_values) == nrow(edges)) edges$label = label_values
    }
    if ("lty" %in% edge_attr_names) {
      line_type = igraph::edge_attr(graph, "lty")
      if (length(line_type) == nrow(edges)) {
        if (is.numeric(line_type)) {
          edges$dashes = line_type != 1
        } else {
          dashed_types = c("dashed", "dotted", "dotdash", "longdash", "twodash")
          edges$dashes = line_type %in% dashed_types
        }
      }
    }
  }
  
  # ============================================================================
  # Create visNetwork
  # ============================================================================
  
  vis_network = visNetwork::visNetwork(nodes, edges)
  vis_network = visNetwork::visOptions(vis_network, highlightNearest = TRUE, nodesIdSelection = TRUE)
  vis_network = visNetwork::visInteraction(vis_network, dragNodes = TRUE, dragView = TRUE, zoomView = TRUE)
  
  # Enable HTML mode if needed
  if (!is.null(label_font)) {
    vis_network = visNetwork::visNodes(vis_network, font = list(multi = "html"))
  }
  
  vis_network
}
