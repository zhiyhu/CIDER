#' @title Initial clustering
#'
#' @description Perform batch-specific initial clustering.
#'
#' @param seu Seurat S4 object. Required.
#' @param batch.var Character. One of the column names of `seu@meta.data`. It
#' is used to partition the Seurat object into smaller ones. Default: "Batch"
#' @param cut.height Numeric. Height used to cut hirerchical trees. Default:
#' 0.4
#' @param nfeatures  Number of high variance genes used. Default: 2000
#' @param additional.vars.to.regress Additional variables to regress out. Needs
#' to among column names of `seu@meta.data`. Default: `NULL`
#' @param dims Number of dimension used for clustering. Passed to Seurat.
#' Default: `1:14`
#' @param resolution Resolution for clustering. Passed to Seurat. Default: 0.6
#' @param verbose Print the progress bar or not. Default: FALSE
#' @param downsampling.size Numeric. The number of cells representing each
#' group. (Default: 40)
#'
#' @return Seurat S4 object with initial cluster information in
#' `initial_cluster` of meta.data.
#' @seealso \code{\link{getIDEr}} \code{\link{finalClustering}}
#' @export
#'
#' @importFrom stats hclust cutree as.dist
#' @import Seurat
#' @importFrom parallel mclapply
#'
initialClustering <- function(seu, batch.var = "Batch",
                              cut.height = 0.4,
                              nfeatures = 2000,
                              additional.vars.to.regress = NULL,
                              dims = seq_len(14), resolution = 0.6,
                              downsampling.size = 50,
                              verbose = FALSE
                              ) {

  if(!batch.var %in% colnames(seu@meta.data)) {
    stop("batch.var does not exist in colnames(seu@meta.data).")
  } else{
    batches <- seu@meta.data[,match(batch.var, colnames(seu@meta.data))]
    if(length(unique(batches)) <= 1){
      stop("Less than 2 batches provided.")
    }
  }

  if(!is.null(additional.vars.to.regress)){
    if(!all(additional.vars.to.regress %in% colnames(seu@meta.data))) {
      stop("additional.vars.to.regress do(es) not
           exist in colnames(seu@meta.data).")
    }
  }

  seu_list <- Seurat::SplitObject(seu, split.by = batch.var)
  seu_list <- mclapply(seu_list, function(x) {
    x <- NormalizeData(x, normalization.method = "LogNormalize",
                       scale.factor = 10000, verbose = verbose)
    x <- FindVariableFeatures(x, selection.method = "vst",
                              nfeatures = nfeatures, verbose = verbose)
    if(is.null(additional.vars.to.regress)){
      x <- ScaleData(x, verbose = verbose)
    } else {
      x <- ScaleData(x, verbose = verbose,
                     vars.to.regress = additional.vars.to.regress)
    }

    x <- RunPCA(x, features = VariableFeatures(object = x), verbose = verbose)
    x <- FindNeighbors(x, dims = dims, verbose = verbose)
    x <- FindClusters(x, resolution = resolution, verbose = verbose)
    return(x)
  })

  dist_coef <- getDistMat(seu_list, downsampling.size = downsampling.size)

  for(seu_itor in seq_len(length(seu_list))){
    tmp <- dist_coef[[seu_itor]] + t(dist_coef[[seu_itor]])
    diag(tmp) <- 1
    tmp <- 1 - tmp
    hc <- hclust(as.dist(tmp), method = "average")
    hres <- cutree(hc, h = cut.height)
    df_hres <- data.frame(hres)
    df_hres$hres <- paste0(df_hres$hres, "_",
                           unique(seu_list[[seu_itor]]$Batch))
    seu_list[[seu_itor]]$inicluster_tmp <-
      paste0(seu_list[[seu_itor]]$seurat_clusters,
             "_", seu_list[[seu_itor]]$Batch)
    seu_list[[seu_itor]]$inicluster <-
      df_hres$hres[match(seu_list[[seu_itor]]$inicluster_tmp,rownames(df_hres))]
  }

  res <- unlist(lapply(seu_list, function(x) return(x$inicluster)))
  res_names <- unlist(lapply(seu_list, function(x) return(colnames(x))))
  seu$initial_cluster <- res[match(colnames(seu), res_names)]

  return(seu)
}

