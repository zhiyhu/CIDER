#' @title Plot distance matrix
#' 
#' 
#' 
#' @export
#' @import pheatmap
plotDistMat <- function(dist_list, use = "coef"){
  
  if(use == "coef"){
    dist_coef <- dist_list[[1]]
  } else if(use == "t"){
    dist_coef <- dist_list[[2]]
  } else if(use == "p"){
    dist_coef <- dist_list[[3]]
  }
  
  
  hm <- list()
  for(i in which(sapply(dist_coef, function(x) return(!is.null(x))))){
    tmp <- dist_coef[[i]] + t(dist_coef[[i]])
    diag(tmp) <- 1
    hm[i] <- pheatmap::pheatmap(tmp, display_numbers = TRUE)
  }
  
  return(hm)
  
}