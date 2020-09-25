#' @param x1 x1
#' @param x2 x2
#' @param method method
#' @importFrom stats cor
measureSimilarity <- function(x1, x2, method = "pearson"){
  if (!is.null(x1) & !is.null(x2)) {
    if(length(x1) != length(x2)) {
      warning("x1 or x2 don't have the same length for similarity measures")
      return(NA)
    }
    if (method %in% c("pearson", "spearman","kendall")) {
      return(cor(x1, x2, method = method))
    } else if (method == "cosine") {
      x <- matrix(cbind(x1, x2))
      return(1-cosine_similarity(x)[1,2])
    }
  } else {
    warning("x1 or x2 are not valid for similarity measures")
    return(NA)
  }
}