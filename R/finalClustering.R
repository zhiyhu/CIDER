
#' @title Final clustering
#' 
#' @description description
#' 
#' @author Zhiyuan Hu
#'   
#' @param seu Seurat S4 object
#' @param dist dist matrix
#' 
#' @return a list
#' 
#' @export
#' 
#' @import stats
#' 
finalClustering <- function(seu, dist, cutree_by = "h", cutree.h = 0.35, cutree.k = 3){
  
  hc <- hclust(as.dist(1-(dist + t(dist)))/2)
  
  if(cutree_by == "h"){
    hcluster <- cutree(hc, h = cutree.h)
  } else {
    hcluster <- cutree(hc, k = cutree.k)
  }
  
  df_merge <-  data.frame(initial_clusters = names(hcluster),
                          final_clusters = hcluster)

  seu$cider.final_cluster <- df_merge$final_cluster[match(seu$cider.inicluster,df_merge$initial_clusters)]
  seu$cider.final_cluster[is.na(seu$cider.final_cluster)] <- seu$inicluster[is.na(seu$cider.final_cluster)]
  
  return(seu)
  
}