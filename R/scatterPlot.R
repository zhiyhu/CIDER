#' @title Plot Network Graph
#' @description Plot network for the initial clusters based on IDEr. The width of edges denotes the similarity between two initial clusters.
#'
#' @param seu Seurat S4 object after the step of `getIDER`. Required.
#' @param reduction list. Output of `getIDER`. Required.
#' @param colour.by A vector of Hex colour codes. If no value is given, a vector of 74 colours will be used. (Default: NULL)
#' @param colvec Character. The type of dist matrix to use. Can be one of "coef", "t" and "p". (Default: coef)
#' @param title Character. It should be one of the colnames of Seurat object meta.data. It is used to colour the vertex of the networkgraph. (Default: NULL)
#' @param sort.by.numbers Numerical. Adjust the thickness of the edges. (Default: 6.5)
#'
#' @details Details
#'
#' @return An igraph object
#'
#' @import ggplot2
#' @importFrom Seurat Reductions
#' @importFrom viridis scale_fill_viridis
#'
#' @export
#'
scatterPlot <- function(seu, reduction, colour.by, colvec = NULL, title = NULL, sort.by.numbers = TRUE, viridis_option = "B") { # function for tSNE plot
  if(is.null(colvec)){
    colvec <- c('#7FC97F','#BEAED4','#FDC086','#FFFF99','#386CB0','#F0027F','#BF5B17','#666666','#1B9E77','#D95F02','#7570B3','#E7298A','#66A61E','#E6AB02','#A6761D','#666666','#A6CEE3','#1F78B4','#B2DF8A','#33A02C','#FB9A99','#E31A1C','#FDBF6F','#FF7F00','#CAB2D6','#6A3D9A','#FFFF99','#B15928','#FBB4AE','#B3CDE3','#CCEBC5','#DECBE4','#FED9A6','#FFFFCC','#E5D8BD','#FDDAEC','#F2F2F2','#B3E2CD','#FDCDAC','#CBD5E8','#F4CAE4','#E6F5C9','#FFF2AE','#F1E2CC','#CCCCCC','#E41A1C','#377EB8','#4DAF4A','#984EA3','#FF7F00','#FFFF33','#A65628','#F781BF','#999999','#66C2A5','#FC8D62','#8DA0CB','#E78AC3','#A6D854','#FFD92F','#E5C494','#B3B3B3','#8DD3C7','#FFFFB3','#BEBADA','#FB8072','#80B1D3','#FDB462','#B3DE69','#FCCDE5','#D9D9D9','#BC80BD','#CCEBC5','#FFED6F')
  }
  
  if(reduction %in% Reductions(seu)){
    df_plot <- data.frame(x = Reductions(seu, reduction)@cell.embeddings[,1],
                          y = Reductions(seu, reduction)@cell.embeddings[,2],
                          stringsAsFactors = FALSE)
  } else {
    stop("Provided reduction name does not exist in the seurat object. Available reduction names can by check by `Reductions(SeuratObject)`.")
  }

  if(!colour.by %in% colnames(seu@meta.data)){
    stop("Provided colour.by name does not exist in the colnames of meta.data of this seurat object.")
  } else {
    idx <- match(colour.by, colnames(seu@meta.data))
    if(mode(seu@meta.data[,idx]) == "numeric" & !is.factor(seu@meta.data[,idx])){ # continues
      df_plot$group <- seu@meta.data[,idx]
    } else {
      df_plot$group <- as.character(seu@meta.data[,idx])
      if(sort.by.numbers) ranked <- names(sort(table(df_plot$group), decreasing = TRUE))
      df_plot$group <- factor(df_plot$group, levels = ranked)
    }
  }
  
  p <- ggplot(df_plot, aes_string(x = "x", y = "y", fill = "group")) + 
    geom_point(alpha = 0.6, size= 0.3)  +  
    theme_classic() + theme(legend.position = "right") + geom_point(col = "grey40", size = 2, shape = 21) +
    theme(axis.ticks.x = element_blank(), axis.ticks.y = element_blank(), axis.text = element_blank())   + 
    xlab(paste0(reduction, "_1")) + ylab(paste0(reduction, "_2")) + labs(fill = colour.by)
  if(!is.null(title)) p <- p + labs(title = title)
  if(mode(seu@meta.data[,idx]) == "numeric" & !is.factor(seu@meta.data[,idx])){
    p <- p + scale_fill_viridis(option = viridis_option) 
  } else {
    p <- p + scale_fill_manual(values = colvec)
  }
  p # plot
}
