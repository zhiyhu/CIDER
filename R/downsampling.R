#' @title Downsampling
#' @description description
#' 
#' @author Zhiyuan Hu
#' 
#' @param metadata metadata
#' @param n.size n.size (default: 50)
#' @param seed seed (default: 12345)
#' @param replace using replace = TRUE if the cluster is smaller than required size
#' @param include using include = TRUE if including the cluster smaller than required size
#' @param lower.cutoff lower cutoff
#' 
#' @return return obj
#' 
#' @export
#' 
#' 
downsampling <- function(metadata, n.size = 50, seed = 12345, include = FALSE, replace = FALSE, lower.cutoff = 3) {
  cluster <- unique(metadata$label)
  tech <- unique(metadata$batch)
  select <- c()
  for(i in cluster){
    for(j in tech) {
      idx <- which(metadata$label %in% i & metadata$batch %in% j)
      if(length(idx) > n.size){
        set.seed(seed)
        select <- c(select, sample(idx, size = n.size, replace = FALSE))
      } else if (length(idx) == n.size) {
        select <- idx
      } else if (length(idx) < n.size & length(idx) >= lower.cutoff) {        
        if(include & !replace) {
          select <- c(select, idx)
        } else if (include & replace){
          set.seed(seed)
          select <- c(select, sample(idx, size = n.size, replace = TRUE))
        }
      }
    }
  }
  return(select)
}
