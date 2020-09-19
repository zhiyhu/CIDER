#' @title Calculate Distance Matrix
#' 
#' @description discription
#' 
#' @author Zhiyuan Hu
#'   
#' @param matrix matrix
#' @param metadata metadata
#' @param verbose print the message and progress bar (default: TRUE)
#' @param method voom or trend
#' @param model one model
#' 
#' @return A list containing
#' 
#' @export
#' 
getDistMat <- function(matrix, metadata, 
                       verbose = TRUE, 
                       method = "voom",
                       model = "one"){
  if(model == "one") {
    dist_list <- calculateDistMatOneModel(matrix = matrix, metadata = metadata, 
                                          verbose = verbose, method = method)
  }
  return(dist_list)
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
#' 
#' @return a list
#' 
#' @export
#' 
#' @import limma edgeR stats
#' 
calculateDistMatOneModel <- function(matrix, metadata, 
                                     verbose = TRUE, 
                                     method = "voom"){
  
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
  
  if("donor" %in% colnames(metadata)){
    df$subb = metadata$donor
    design <- model.matrix(~  0 + g + subb + detrate, data = df)
  } else {
    design <- model.matrix(~  0 + g + detrate, data = df)
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
    fit <- lmFit(v, design)
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
