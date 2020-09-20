#' @title Downsampling cells
#' 
#' @description Downsampling cells from each group 
#' 
#' @param metadata Data frame. It includes at least 2 columns, label and batch. Each row corresponds to one cell. Required.
#' @param n.size Numeric. The number of cells used in each group. (Default: 35)
#' @param seed Numeric. Seed used to sample. (Default: 12345)
#' @param include Boolean. Using `include = TRUE` to include the group smaller than required size. (Default: FALSE)
#' #' @param replace Boolean. Using `replace = TRUE` if the group is smaller than required size and some cells will be repeatedly used. (Default: FALSE)
#' @param lower.cutoff Numeric. The minimum size of groups to keep. (Default: 3)
#' 
#' @return A numeric list of which cells will be kept for downstream computation.
#' 
#' @export
#' 
downsampling <- function(metadata, n.size = 35, seed = 12345, include = FALSE, replace = FALSE, lower.cutoff = 3) {
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
