#' @title Plot Similarity Matrix
#'
#' @description description
#'
#' @param dist.list List of list of similarity matrix. Output of function `getDistMat()`. Required.
#' @param use Choose from "coef", "t" and "p". (Default: coef)
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
  } else if (use == "t") {
    dist_coef <- dist.list[[2]]
  } else if (use == "p") {
    dist_coef <- dist.list[[3]]
  }

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
