
#' @title Final clustering
#' 
#' @description description
#' 
#' @author Zhiyuan Hu
#'   
#' @param seu Seurat S4 object
#' @param dist dist matrix
#' @param cutree.by cutree.by
#' @param cutree.h cutree.h
#' @param cutree.k cutree.k
#' @param hc.method hc.method
#' 
#' @return a list
#' 
#' @export
#' 
#' @import stats
#' 
finalClustering <- function(seu, dist, cutree.by = "h", 
                            cutree.h = 0.35, cutree.k = 3,
                            hc.method = "complete"){
  
  hc <- hclust(as.dist(1-(dist + t(dist)))/2, method = hc.method)
  
  if(cutree.by == "h"){
    hcluster <- cutree(hc, h = cutree.h)
  } else {
    hcluster <- cutree(hc, k = cutree.k)
  }
  
  df_merge <-  data.frame(initial_clusters = names(hcluster),
                          final_clusters = hcluster)
  
  seu$cider_final_cluster <- df_merge$final_cluster[match(seu$initial_cluster,
                                                          df_merge$initial_clusters)]
  seu$cider_final_cluster[is.na(seu$cider_final_cluster)] <- seu$initial_cluster[is.na(seu$cider_final_cluster)]
  
  return(seu)
  
}