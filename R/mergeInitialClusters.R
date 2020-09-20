#' @title Merge Initial Clusters
#' 
#' @param seu_list seu_list
#' @param dist_list dist_list
#' @param use use
#' @param method method
#' @param hc.method hc.method
#' @param cutree.by cutree.by
#' @param cutree.h cutree.h
#' @param cutree.k cutree.k 
#' 
#' @return a list of seurat objects
#' 
#' @export
#' 
#' @importFrom stats cutree hclust
#' @import Seurat
#' 
#' 
#' 
mergeInitialClusters <- function(seu_list, dist_list, use = "coef", method = "hc", 
                                 hc.method = "average", cutree.by = "h", cutree.h = 0.6, cutree.k = 3){
  
  if(use == "coef"){
    dist_coef <- dist_list[[1]]
  } else if(use == "t"){
    dist_coef <- dist_list[[2]]
  } else if(use == "p"){
    dist_coef <- dist_list[[3]]
  }
  
  for(seu_itor in 1:length(seu_list)){
    
    hc <- hclust(1-as.dist(dist_coef[[seu_itor]] + t(dist_coef[[seu_itor]])), method = hc.method)
    
    if(cutree.by == "h"){
      hres <- cutree(hc, h = cutree.h)
    } else {
      hres <- cutree(hc, k = cutree.k)
    }
    
    df_hres <- data.frame(hres)
    df_hres$hres <- paste0(df_hres$hres, "_", unique(seu_list[[seu_itor]]$Batch))
    seu_list[[seu_itor]]$inicluster_tmp <- paste0(seu_list[[seu_itor]]$seurat_clusters, "_", seu_list[[seu_itor]]$Batch)
    seu_list[[seu_itor]]$inicluster <- df_hres$hres[match(seu_list[[seu_itor]]$inicluster_tmp,rownames(df_hres))]
  }
  return(seu_list)
}

#' @title Gather initial cluster names
#' 
#' @param seu_list seu_lisst
#' @param seu seu
#' 
#' @return a seurat object
#' 
#' @import Seurat 
#' @export
gatherInitialClusters <- function(seu_list, seu){
  
  tmp <- unlist(sapply(seu_list, function(x) return(x$inicluster_tmp)))
  names(tmp) <- unlist(sapply(seu_list, function(x) return(colnames(x@assays$RNA@counts))))
  seu$initial_cluster_tmp <- tmp[match(colnames(seu), names(tmp))]
  
  tmp <- unlist(sapply(seu_list, function(x) return(x$inicluster)))
  names(tmp) <- unlist(sapply(seu_list, function(x) return(colnames(x@assays$RNA@counts))))
  seu$initial_cluster <- tmp[match(colnames(seu), names(tmp))]
  
  return(seu)
}
