#' @title Calculate Distance Matrix
#' 
#' @description discription
#' 
#' @author Zhiyuan Hu
#'   
#' @param seu_list seu list
#' @param verbose print the message and progress bar (default: TRUE)
#' @param method voom or trend
#' @param model one model
#' @param additional.variate additional.variate
#' @param downsampling.size downsampling.size
#' @param downsampling.include downsampling.include
#' @param downsampling.replace downsampling.replace
#' 
#' @return A list containing
#' 
#' @export
#' 
#' @import Seurat utils
#' 
getDistMat <- function(seu_list, 
                       verbose = TRUE, 
                       method = "trend",
                       model = "one",
                       additional.variate = "donor",
                       downsampling.size = 35,
                       downsampling.include = TRUE,
                       downsampling.replace = TRUE){

  dist_p <- list()
  dist_t <- list()
  dist_coef <- list()
  
  # create progress bar
  if(verbose == TRUE) {
    pb <- txtProgressBar(min = 0, max = length(seu_list), style = 3)
    k <- 1
  }
  
  for(seu_itor in 1:length(seu_list)){
    
    df_info <- data.frame(label = seu_list[[seu_itor]]$seurat_clusters,
                          batch = seu_list[[seu_itor]]$Batch,
                          donor = seu_list[[seu_itor]]$Tissue)
    
    idx <- downsampling(metadata = df_info, n.size = downsampling.size, 
                        include = downsampling.include, replace = downsampling.replace)
    idx <- sort(idx)
    
    to_add <-  idx[duplicated(idx)]
    idx <- idx[!duplicated(idx)] 
    matrix <- as.matrix(seu_list[[seu_itor]]@assays$RNA@counts[,idx])
    
    if(length(to_add) > 0) {
      
      matrix2 <- data.frame(seu_list[[seu_itor]]@assays$RNA@counts[,to_add])
      colnames(matrix2) <- paste0(colnames(matrix2), 1:ncol(matrix2))
      matrix2 <- as.matrix(matrix2)
      matrix <- cbind(matrix, matrix2)
      rm(matrix2)
      
    }
    
    if(model == "one") {
      
      if(length(unique(df_info$label[idx])) > 2)
        dist_list <- calculateDistMatOneModel(matrix = matrix, metadata = df_info[c(idx, to_add),], 
                                              verbose = verbose, method = method,
                                              additional.variate = additional.variate)
    }
    
    dist_coef[[seu_itor]] <- dist_list[[1]]
    dist_t[[seu_itor]] <- dist_list[[2]]
    dist_p[[seu_itor]] <- dist_list[[3]]
    
    if(verbose == TRUE) {
      setTxtProgressBar(pb, k) # progress bar
      k <- k+1
    }
    
  }
  
  if(verbose == TRUE) {
    close(pb) # close progress bar
  }
  
  return(list(dist_coef,dist_t, dist_p))
}



#' @title Calculate distance matrix with in one model
#' 
#' @description description
#' 
#' @author Zhiyuan Hu
#'   
#' @param matrix matrix
#' @param metadata metadata
#' @param verbose print the message and progress bar (default: TRUE)
#' @param method voom or trend
#' @param additional.variate additional.variate
#' 
#' @return a list
#' 
#' @export
#' 
#' @import limma edgeR stats
#' 
calculateDistMatOneModel <- function(matrix, metadata, 
                                     verbose = TRUE, 
                                     method = "voom",
                                     additional.variate = "donor"){
  
  keep <- rowSums(matrix > 0.5) > 5 
  dge <- DGEList(counts = matrix[keep,,drop=F]) # make a edgeR object
  dge <- dge[!grepl("ERCC-", rownames(dge)),] # remove ERCC
  dge <- dge[!grepl("MT-", rownames(dge)),]
  
  df <- data.frame(g = paste(metadata$label, metadata$batch, sep = "_"),
                   b = metadata$batch, ## batch
                   c = metadata$label, stringsAsFactors = F) ## label
  
  df$detrate <- scale(colMeans(matrix > 0))[,1]
  rownames(df) <- colnames(matrix)
  
  N <- length(unique(df$g))
  
  combinations <- data.frame(g1 = rep(unique(df$g), each = N), g2 = rep(unique(df$g), N), stringsAsFactors = FALSE)
  combinations <- combinations[combinations$g1 != combinations$g2, ]
  combinations <- sort(combinations)
  
  idx <- c()
  for(i in 2:nrow(combinations)){
    if(!combinations$g2[i] %in% combinations$g1[1:(i-1)]) {
      idx <- c(idx, i)
    }
  }
  
  combinations <- combinations[c(1,idx),]
  rownames(combinations) <- 1:nrow(combinations)
  
  dist_p <- dist_t <- dist_coef <- matrix(0, nrow = N, ncol = N)
  colnames(dist_p) <- rownames(dist_p) <- sort(unique(df$g))
  colnames(dist_t) <- rownames(dist_t) <- sort(unique(df$g))
  colnames(dist_coef) <- rownames(dist_coef) <- sort(unique(df$g))
  
  
  if(!is.null(additional.variate) & additional.variate %in% colnames(metadata)){
    df$subb <- metadata[,colnames(metadata) == additional.variate]
    design <- model.matrix(~  0 + g + subb + detrate, data = df)
    message("model used: ~  0 + g + subb + detrate")
  } else {
    design <- model.matrix(~  0 + g + detrate, data = df)
    message("model used: ~  0 + g + detrate")
  }
  
  groups <- sort(unique(paste0("g", df$g)))
  n_groups <- length(groups)
  
  df_contrasts <- data.frame(target_group = groups, contrast = NA)
  
  for(i in 1:n_groups){
    df_contrasts$contrast[i] <- paste0(groups[i], "-(", paste(groups[-i], collapse = "+"), ")/", (n_groups-1))
  }
  # contrast matrix
  contrast_m <- makeContrasts(contrasts = df_contrasts$contrast, levels = design)
  colnames(contrast_m) <- groups
  
  if (method == "voom") {
    v <- voom(dge, design, plot = FALSE)
    suppressMessages(fit <- lmFit(v, design))
    group_fit <- contrasts.fit(fit, contrast_m) 
    group_fit <- eBayes(group_fit)
    
  } else if (method == "trend") {
    logCPM <- cpm(dge, log=TRUE, prior.count=3)
    fit <- lmFit(logCPM, design)
    fit <- contrasts.fit(fit, contrast_m)
    group_fit <- eBayes(fit, trend=TRUE, robust = TRUE)
    
  }
  
  for(i in 1:ncol(group_fit$p.value)) {
    group_fit$p.value[,i] <- group_fit$p.value[,i] + 0.00000001
  }
  
  # pairwise comparison
  for(i in 1:nrow(combinations)){
    
    idx1 <- rownames(dist_coef) == combinations$g1[i]
    idx2 <- colnames(dist_coef) == combinations$g2[i]
    
    pos1 <- df_contrasts$target_group == paste0("g", combinations$g1[i])
    pos2 <- df_contrasts$target_group == paste0("g", combinations$g2[i])
    
    dist_coef[idx1, idx2] <- cor(coef(group_fit)[,pos1], coef(group_fit)[,pos2])
    dist_t[idx1, idx2]    <- cor(group_fit$t[,pos1], group_fit$t[,pos2])
    dist_p[idx1, idx2]    <- cor(-log10(group_fit$p.value)[,pos1]*sign(coef(group_fit)[,pos1]), 
                                 -log10(group_fit$p.value)[,pos2]*sign(coef(group_fit)[,pos2]), 
                                 use = "pairwise.complete.obs")
    
  }
  
  return(list(dist_coef, dist_t, dist_p))
}
