#' @title Plot Network Graph
#' @description Plot network for the initial clusters based on IDEr. The width of edges denotes the similarity between two initial clusters.
#'
#' @param seu Seurat S4 object after the step of `getIDER`. Required.
#' @param dist A list. Output of `getIDER`. Required.
#' @param col.vector A vector of Hex colour codes. If no value is given, a vector of 74 colours will be used. (Default: NULL)
#' @param use Character. The type of dist matrix to use. Can be one of "coef", "t" and "p". (Default: coef)
#' @param colour.by Character. It should be one of the colnames of Seurat object meta.data. It is used to colour the vertex of the networkgraph. (Default: NULL)
#' @param weight.factor Numerical. Adjust the thickness of the edges. (Default: 6.5)
#' @param vertex.size Numerical. Adjsut the size of vertexes. (Default: 1)
#'
#' @details Details
#'
#' @return An igraph object
#'
#' @seealso \code{\link{getIDEr}} \code{\link[igraph]{graph_from_data_frame}}
#'
#' @import Seurat igraph
#' @importFrom graphics plot legend
#'
#'
#' @export
#'
plotNetwork <- function(seu, dist, col.vector = NULL, use = "coef",
                        colour.by = NULL, weight.factor = 6.5, vertex.size = 1) {
  select <- dist[[2]]

  if (!is.null(colour.by)) {
    if (colour.by %in% colnames(seu@meta.data)) {
      seu$Group <- seu@meta.data[, colnames(seu@meta.data) == colour.by]
    } else {
      warning("'colour.by' is not in the colnames of Seurat object meta.data. Please check.")
    }
  } else {
    seu$Group <- 1
  }

  df <- data.frame(
    g = seu$initial_cluster[select],
    b = seu$Batch[select], ## batch
    stringsAsFactors = F
  ) ## label

  N <- length(unique(df$g)) # number of groups

  # get the dataframe of combinations/pairs for comparison
  combinations <- dist[[3]]

  df_plot <- data.frame(
    g = seu$Group[select], # colour.by
    b = seu$Batch[select], # batch
    c = seu$initial_cluster[select], # label; initial cluster
    stringsAsFactors = F
  )

  df_plot$combination <- paste0(df_plot$g, "-", df_plot$c)
  freq <- table(df_plot$combination)
  df_plot$freq <- freq[match(df_plot$combination, names(freq))]
  df_plot <- df_plot[order(df_plot$freq, decreasing = TRUE), ]
  df_plot <- unique(df_plot)

  edges <- data.frame(from = combinations$g1, to = combinations$g2, weight = NA) # edges

  if (use == "coef") { # type of dist matrix to use
    tmp <- dist[[1]] + t(dist[[1]])
  } else if (use == "t") {
    tmp <- dist[[2]] + t(dist[[2]])
  } else if (use == "p") {
    tmp <- dist[[3]] + t(dist[[3]])
  } else {
    warning("'use' is incorrect. Please use one of coef, t or p.")
  }

  for (i in 1:nrow(edges)) {
    edges$weight[i] <- tmp[rownames(tmp) == edges$from[i], colnames(tmp) == edges$to[i]]
  }

  edges$weight[edges$weight < 0] <- 0
  edges <- edges[edges$weight > 0, ]
  net <- igraph::graph_from_data_frame(edges, directed = FALSE)
  E(net)$width <- 2^(E(net)$weight * weight.factor)
  vg_names <- attr(V(net), "names")

  if (is.null(col.vector)) {
    col.vector <- c(
      "#7FC97F", "#BEAED4", "#FDC086", "#FFFF99", "#386CB0", "#F0027F", "#BF5B17", "#666666", "#1B9E77",
      "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02", "#A6761D", "#666666", "#A6CEE3", "#1F78B4",
      "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C", "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A", "#FFFF99",
      "#B15928", "#FBB4AE", "#B3CDE3", "#CCEBC5", "#DECBE4", "#FED9A6", "#FFFFCC", "#E5D8BD", "#FDDAEC",
      "#F2F2F2", "#B3E2CD", "#FDCDAC", "#CBD5E8", "#F4CAE4", "#E6F5C9", "#FFF2AE", "#F1E2CC", "#CCCCCC",
      "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33", "#A65628", "#F781BF", "#999999",
      "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", "#8DD3C7",
      "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5", "#D9D9D9", "#BC80BD",
      "#CCEBC5", "#FFED6F"
    )
  }

  if (length(unique(df_plot$g)) <= length(col.vector) & length(unique(df_plot$g)) > 1) {
    df_plot_cols <- data.frame(
      group = unique(df_plot$g),
      col = col.vector[1:length(unique(df_plot$g))], stringsAsFactors = FALSE
    )
    df_plot$col <- df_plot_cols$col[match(df_plot$g, df_plot_cols$group)]
    V(net)$color <- df_plot$col[match(vg_names, df_plot$c)]
  } else if (length(unique(df_plot$g)) > length(col.vector)) {
    warning("The number of provided colours is less than needed. So no colour is used.")
  }

  V(net)$frame.color <- "#777777"
  V(net)$size <- 20 * vertex.size
  V(net)$label.family <- "Helvetica"

  if (length(unique(df_plot$g)) > 1) {
    plot(net); legend(
      x = -1.5, y = -1.1, df_plot_cols$group, pch = 21,
      col = "#777777", pt.bg = df_plot_cols$col, pt.cex = 2, cex = .8, bty = "n", ncol = 2
    )
  } else {
    plot(net)
  }
  return(net)
}
