#' @title Merge Initial Clusters
#'
#' @param seu_list A list containing Seurat objects. Required.
#' @param dist_list A list containing similarity matrices. The output of
#' `getDistMat ()`
#' @param use Default: "coef". No other option available currently.
#' @param method method = "hc"
#' @param hc.method Passed to the `method` parameter of `hclust()`. Default:
#' "average"
#' @param cutree.by Cut trees by height ("h", default) or number of
#' clusters ("k")
#' @param cutree.h Height used to cut the tree. Default: 0.6.
#' @param cutree.k Number of clusters used to cut the tree. Default: 3.
#'
#' @return a list of Seurat objects containing the updated initial clustering
#' information in `seu_list[[seu_itor]]$inicluster`. The original initial
#' cluster information is stored in `seu_list[[seu_itor]]$inicluster_tmp`.
#' @seealso \code{\link{hclust}}  \code{\link{cutree}}
#' \code{\link{gatherInitialClusters}} \code{\link{initialClustering}}
#' @export
#'
#' @importFrom stats cutree hclust as.dist
#' @import Seurat
#'
#'
mergeInitialClusters <- function(seu_list, dist_list, use = "coef",
                                 method = "hc",
                                 hc.method = "average", cutree.by = "h",
                                 cutree.h = 0.6, cutree.k = 3) {
  if (use == "coef") {
    dist_coef <- dist_list[[1]]
  } else if (use == "t") {
    dist_coef <- dist_list[[2]]
  } else if (use == "p") {
    dist_coef <- dist_list[[3]]
  }

  for (seu_itor in seq_len(length(seu_list))) {
    hc <- hclust(1 - as.dist(dist_coef[[seu_itor]] + t(dist_coef[[seu_itor]])),
                 method = hc.method)

    if (cutree.by == "h") {
      hres <- cutree(hc, h = cutree.h)
    } else {
      hres <- cutree(hc, k = cutree.k)
    }

    df_hres <- data.frame(hres)
    df_hres$hres <- paste0(df_hres$hres, "_",
                           unique(seu_list[[seu_itor]]$Batch))
    seu_list[[seu_itor]]$inicluster_tmp <-
      paste0(seu_list[[seu_itor]]$seurat_clusters, "_",
             seu_list[[seu_itor]]$Batch)
    seu_list[[seu_itor]]$inicluster <-
      df_hres$hres[match(seu_list[[seu_itor]]$inicluster_tmp,
                         rownames(df_hres))]
  }
  return(seu_list)
}

#' @title Gather initial cluster names
#' @describeIn Merge initial clustering results from a Seurat object list to one
#' Seurat object. Follows the function `mergeInitialClusters`.
#' @param seu_list  A list containing Seurat objects. Required.
#' @param seu A Seurat object
#'
#' @return A Seurat object containing initial clustering  results in
#' `seu$initial_cluster`.
#' @seealso \code{\link{mergeInitialClusters}}
#' @import Seurat
#' @export
gatherInitialClusters <- function(seu_list, seu) {
  tmp <- unlist(vapply(seu_list, function(x) {
    return(x$inicluster_tmp)
  }))
  names(tmp) <- unlist(vapply(seu_list, function(x) {
    return(colnames(x@assays$RNA@counts))
  }))
  seu$initial_cluster_tmp <- tmp[match(colnames(seu), names(tmp))]

  tmp <- unlist(vapply(seu_list, function(x) {
    return(x$inicluster)
  }))
  names(tmp) <- unlist(vapply(seu_list, function(x) {
    return(colnames(x@assays$RNA@counts))
  }))
  seu$initial_cluster <- tmp[match(colnames(seu), names(tmp))]

  return(seu)
}
