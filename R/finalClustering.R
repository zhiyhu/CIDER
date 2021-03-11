#' @title Final clustering
#'
#' @description Merge initial clusters into final clusters based on the matrix of IDEr.
#'
#' @param seu Seurat S4 object after the step of `getIDEr`. Required.
#' @param dist A list. Output of `getIDEr`. Required.
#' @param use Character. The type of dist matrix to use. Can be one of "coef" and "p". (Default: coef)
#' @param cutree.by Character. Cut the tree by which parameter, height ("h") or number of clusters ("k"). (Default: h)
#' @param cutree.h Numeric between 0 and 1. The height used to cut the tree. Ignored if `cutree.by = 'k`. (Default: 0.45)
#' @param cutree.k Numeric/integer. Used to cut the tree. Ignored if `cutree.by = 'h`. (Default: 3)
#' @param hc.method Character. Used to choose the hierarchical clustering method.
#'
#' @return Seurat S4 object with final clustering results in `cider_final_cluster` of meta.data.
#'
#' @seealso \code{\link{getIDEr}} \code{\link[stats]{hclust}}
#'
#' @export
#'
#' @importFrom stats hclust cutree as.dist
#'
finalClustering <- function(seu, dist, use = "coef",
                            cutree.by = "h", cutree.h = 0.45, cutree.k = 3,
                            hc.method = "complete") {
  if (use == "coef") {
    tmp <- dist[[1]] + t(dist[[1]])
  } else if (use == "p") {
    tmp <- dist[[2]] + t(dist[[2]])
  } else {
    warning("'use' is incorrect. Please use one of coef or p.")
  }

  hc <- hclust(as.dist(1 - tmp) / 2, method = hc.method)

  if (cutree.by == "h") {
    hcluster <- cutree(hc, h = cutree.h)
  } else {
    hcluster <- cutree(hc, k = cutree.k)
  }

  df_merge <- data.frame(
    initial_clusters = names(hcluster),
    final_clusters = hcluster
  )

  seu$cider_final_cluster <- df_merge$final_cluster[match(
    seu$initial_cluster,
    df_merge$initial_clusters
  )]
  seu$cider_final_cluster[is.na(seu$cider_final_cluster)] <- seu$initial_cluster[is.na(seu$cider_final_cluster)]

  return(seu)
}
