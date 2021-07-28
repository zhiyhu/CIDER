#' @title Plot Similarity Matrix
#'
#' @description description
#'
#' @param dist.list List of list of similarity matrix. Output of function `getDistMat()`. Required.
#' @param use Choose from "coef" and "p". (Default: coef)
#'
#' @export
#'
#' @import pheatmap
#'
#' @seealso \code{\link{getDistMat}}
#'
plotDistMat <- function(dist.list, use = "coef") {
  if (use == "coef") {
    dist_coef <- dist.list[[1]]
  } else if (use == "p") {
    dist_coef <- dist.list[[2]]
  }
  
  # else if (use == "t") {
  #   dist_coef <- dist.list[[2]]

  hm <- list()
  for (i in which(sapply(dist_coef, function(x) {
    return(!is.null(x))
  }))) {
    tmp <- dist_coef[[i]] + t(dist_coef[[i]])
    diag(tmp) <- 1
    hm[i] <- pheatmap::pheatmap(tmp, display_numbers = TRUE)
  }
  return(hm)
}


#' @import pheatmap
#' @export
plotHeatmap <- function(seu, ider) {
  idx <- getSharedGroups(seu, ider[[1]])
  shared_g <- idx[[1]]
  idx1 <- idx[[2]]
  idx2 <- idx[[3]]
  
  pheatmap::pheatmap(
    ider[[1]][idx1, idx2],
    border_color = "grey20",
    color = viridis::inferno(50),
    display_numbers = TRUE,
    cluster_rows = FALSE,
    cluster_cols = FALSE
  )
  
}