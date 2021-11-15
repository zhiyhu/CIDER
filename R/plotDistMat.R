#' @title Plot Similarity Matrix with pheatmap
#'
#' @param dist.list Output of function `getDistMat()`. Required.
#' @param use Default: "coef". No other option currently that can be used.
#' @return A pheatmap showing the similarity matrix
#' @export
#'
#' @import pheatmap
#' @seealso \code{\link{getDistMat}}
plotDistMat <- function(dist.list, use = "coef") {
  if (use == "coef") {
    dist_coef <- dist.list[[1]]
  } else if (use == "p") {
    dist_coef <- dist.list[[2]]
  }

  hm <- list()
  for (i in which(vapply(dist_coef, function(x) {
    return(!is.null(x))
  }))) {
    tmp <- dist_coef[[i]] + t(dist_coef[[i]])
    diag(tmp) <- 1
    hm[i] <- pheatmap::pheatmap(tmp, display_numbers = TRUE)
  }
  return(hm)
}

#' @title Plot Heatmap for the IDER-based similarity matrix
#'
#' @param seu An Seurat object.
#' @param ider Output of function `getIDEr`.
#' @return A heatmap shows the similarity between shared groups in two batches
#'
#' @seealso \code{\link{getIDEr}}
#' @import pheatmap viridis
#' @export
#' @examples
#' \dontrun{
#'   plotHeatmap(seu, ider)
#' }
plotHeatmap <- function(seu, ider) {
  idx <- getSharedGroups(seu, ider[[1]])
  shared_g <- idx[[1]] # shared groups
  if(length(shared_g) == 1){
    stop("No shared groups.")
  }
  idx1 <- idx[[2]] # rownames
  idx2 <- idx[[3]] # colnames

  pheatmap::pheatmap(
    ider[[1]][idx1, idx2],
    border_color = "grey20",
    color = viridis::inferno(50),
    display_numbers = TRUE,
    cluster_rows = FALSE,
    cluster_cols = FALSE
  )
}
