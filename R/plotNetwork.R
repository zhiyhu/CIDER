#' @title plot Network
#' @description plot network
#'
#' @param seu seu
#' @param dist dist
#' @param col_vector col_vector
#' @param use use
#' 
#' @import Seurat igraph
#' 
#' 
plotNetwork <- function(seu, dist, col_vector, use = "coef"){
  
  select <- dist[[4]]
  
  df <- data.frame(g = seu$initial_cluster[select],
                   b = seu$Batch[select], ## batch
                   ground_truth = seu$Group[select],
                   stringsAsFactors = F) ## label
  
  N <- length(unique(df$g)) # number of groups
  
  # get the dataframe of combinations/pairs for comparison
  combinations <- data.frame(g1 = rep(unique(df$g), each = N), g2 = rep(unique(df$g), N), stringsAsFactors = FALSE)
  combinations <- combinations[combinations$g1 != combinations$g2, ]
  combinations$b1 <- df$b[match(combinations$g1, df$g)]
  combinations$b2 <- df$b[match(combinations$g2, df$g)]
  combinations <- combinations[combinations$b1!=combinations$b2,]
  idx <- c()
  for(i in 2:nrow(combinations)){
    if(!combinations$g2[i] %in% combinations$g1[1:(i-1)]) {
      idx <- c(idx, i)
    }
  }
  combinations <- combinations[c(1,idx),]
  rownames(combinations) <- 1:nrow(combinations)
  
  df_plot <- data.frame(g = seu$Group[select],
                        b = seu$Batch[select], ## batch
                        c = seu$initial_cluster[select], ## label
                        stringsAsFactors = F) 
  
  df_plot$combination <- paste0(df_plot$g, "-", df_plot$c)
  freq <- table(df_plot$combination)
  df_plot$freq <- freq[match(df_plot$combination, names(freq))]
  df_plot <- df_plot[order(df_plot$freq, decreasing = TRUE),]
  df_plot <- unique(df_plot)
  
  N <- length(unique(df_plot$g))
  
  edges <- data.frame(from = combinations$g1, to = combinations$g2, weight = NA) # edges
  
  if(use == "coef"){
    tmp <- dist[[1]] + t(dist[[1]])
  }
  
  for(i in 1:nrow(edges)) {
    edges$weight[i] <- tmp[rownames(tmp) == edges$from[i], colnames(tmp) == edges$to[i]]
  }
  
  edges$weight[edges$weight < 0] <- 0
  edges <- igraph::edges[edges$weight > 0, ]
  net <- igraph::graph_from_data_frame(edges, directed = FALSE)
  E(net)$width <- 2^(E(net)$weight * 6)
  vg_names <- attr(V(net), "names")
  df_plot_cols <- data.frame(group = unique(df_plot$g),
                             col = col_vector[1:length(unique(df_plot$g))], stringsAsFactors = FALSE)
  df_plot$col <- df_plot_cols$col[match(df_plot$g, df_plot_cols$group)]
  
  V(net)$color <- df_plot$col[match(vg_names, df_plot$c)]
  V(net)$frame.color <- "#777777"
  V(net)$size <- 10
  V(net)$label.family <- "Helvetica"
  
  p <- plot(net);legend(x=-1.5, y=-1.1, df_plot_cols$group, pch=21,
                   col="#777777", pt.bg=df_plot_cols$col, pt.cex=2, cex=.8, bty="n", ncol=1)
  
  return(p)
}