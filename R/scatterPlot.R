#' @title Scatterplot by a selected feature
#' @description Scatterplot of a Seurat object based on dimension reduction.
#'
#' @param seu Seurat S4 object after the step of `getIDER`. Required.
#' @param reduction Character. The dimension reduction used to plot. Common
#' options: "pca", "tsne", "umap". The availability of dimension reduction
#' can be checked by `Reductions(seu)`.
#' @param colour.by Character. One of the column names of `seu@meta.data`.
#' Can be either discreet or continuous variables.
#' @param colvec A vector of Hex colour codes. If no value is given (default),
#' a vector of 74 colours will be used.
#' @param title Character. Title of the figure.
#' @param sort.by.numbers Boolean. Whether to sort the groups by the number
#'  of cells.(Default: True)
#' @param viridis_option viridis_option. (Default: B)
#' @return a scatter plot
#' @import ggplot2
#' @importFrom Seurat Reductions
#' @importFrom viridis scale_fill_viridis
#' @export
#'
#' @examples
#' \dontrun{
#' library(Seurat)
#' seu <- NormalizeData(seu, verbose = FALSE) # input Seurat object
#' seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 2000,
#'  verbose = FALSE)
#' seu <- ScaleData(seu, verbose = FALSE)
#' seu <- RunPCA(seu, npcs = 20, verbose = FALSE)
#' scatterPlot(seu, "pca",colour.by = "Batch", title = "PCA")
#' }
scatterPlot <- function(seu, reduction, colour.by, colvec = NULL,
                        title = NULL, sort.by.numbers = TRUE,
                        viridis_option = "B") {
  # function for tSNE plot
  if(is.null(colvec)){
    colvec <-
      c('#7FC97F','#BEAED4','#FDC086','#FFFF99','#386CB0',
        '#F0027F','#BF5B17','#666666','#1B9E77','#D95F02',
        '#7570B3','#E7298A','#66A61E','#E6AB02','#A6761D',
        '#666666','#A6CEE3','#1F78B4','#B2DF8A','#33A02C',
        '#FB9A99','#E31A1C','#FDBF6F','#FF7F00','#CAB2D6',
        '#6A3D9A','#FFFF99','#B15928','#FBB4AE','#B3CDE3',
        '#CCEBC5','#DECBE4','#FED9A6','#FFFFCC','#E5D8BD',
        '#FDDAEC','#F2F2F2','#B3E2CD','#FDCDAC','#CBD5E8',
        '#F4CAE4','#E6F5C9','#FFF2AE','#F1E2CC','#CCCCCC',
        '#E41A1C','#377EB8','#4DAF4A','#984EA3','#FF7F00',
        '#FFFF33','#A65628','#F781BF','#999999','#66C2A5',
        '#FC8D62','#8DA0CB','#E78AC3','#A6D854','#FFD92F',
        '#E5C494','#B3B3B3','#8DD3C7','#FFFFB3','#BEBADA',
        '#FB8072','#80B1D3','#FDB462','#B3DE69','#FCCDE5',
        '#D9D9D9','#BC80BD','#CCEBC5','#FFED6F')
  }

  if(reduction %in% Reductions(seu)){
    df_plot <- data.frame(x = Reductions(seu, reduction)@cell.embeddings[,1],
                          y = Reductions(seu, reduction)@cell.embeddings[,2],
                          stringsAsFactors = FALSE)
  } else {
    stop("Provided reduction name does not exist in the seurat object.
         Available reduction names can by check by `Reductions(SeuratObject)`.")
  }

  if(!colour.by %in% colnames(seu@meta.data)){
    stop("Provided colour.by name does not exist in
         the colnames of meta.data of this seurat object.")
  } else {
    idx <- match(colour.by, colnames(seu@meta.data))
    if(mode(seu@meta.data[,idx]) == "numeric" &
       !is.factor(seu@meta.data[,idx])){ # continues
      df_plot$group <- seu@meta.data[,idx]
    } else {
      df_plot$group <- as.character(seu@meta.data[,idx])
      if(sort.by.numbers) ranked <- names(sort(table(df_plot$group),
                                               decreasing = TRUE))
      df_plot$group <- factor(df_plot$group, levels = ranked)
    }
  }

  p <- ggplot(df_plot, aes_string(x = "x", y = "y", fill = "group")) +
    geom_point(alpha = 0.3, size= 0.1)  +
    theme_classic() + theme(legend.position = "right") +
    geom_point(col = "grey40", size = 2, shape = 21) +
    theme(axis.ticks.x = element_blank(), axis.ticks.y = element_blank(),
          axis.text = element_blank())   +
    xlab(paste0(reduction, "_1")) + ylab(paste0(reduction, "_2")) +
    labs(fill = colour.by)
  if(!is.null(title)) p <- p + labs(title = title)
  if(mode(seu@meta.data[,idx]) == "numeric" & !is.factor(seu@meta.data[,idx])){
    p <- p + scale_fill_viridis(option = viridis_option)
  } else {
    p <- p + scale_fill_manual(values = colvec)
  }
  p # plot
}
