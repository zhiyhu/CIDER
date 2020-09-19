#' @title Get IDEr matrix
#' 
#' @description description
#' 
#' @author Zhiyuan Hu
#'   
#' @param seu Seurat S4 object
#' @param verbose print the message and progress bar (default: TRUE)
#' @param method voom or trend
#' 
#' @return a list
#' 
#' @export
#' 
#' @import limma edgeR stats
#' 
getIDEr <- function(seu, 
                    downsampling.size = 40,
                    downsampling.include = TRUE,
                    downsampling.replace = TRUE,
                    verbose = TRUE, 
                    bg.downsampling.factor = 1, 
                    method = "voom"){
  
  
  ## merge seurat list
  metadata <- data.frame(label = seu$initial_cluster,
                         batch = seu$Batch,
                         ground_truth = seu$Group, stringsAsFactors = FALSE)
  
  select <- downsampling(metadata = metadata, n.size = downsampling.size, 
                         include = downsampling.include, replace = downsampling.replace)
  
  matrix <- as.matrix(seu@assays$RNA@counts[,select])
  colnames(matrix) <- paste0(colnames(matrix),1:ncol(matrix)) # avoid duplication
  keep <- rowSums(matrix > 0.5) > 5 
  dge <- edgeR::DGEList(counts = matrix[keep,,drop=F]) # make a edgeR object
  dge <- dge[!grepl("ERCC-", rownames(dge)),] # remove ERCC
  dge <- dge[!grepl("MT-", rownames(dge)),]

  df <- data.frame(g = metadata$label[select],
                   b = metadata$batch[select], ## batch
                   ground_truth = metadata$ground_truth[select],
                   stringsAsFactors = F) ## label
  
  df$detrate <- scale(colMeans(matrix > 0))[,1]
  rownames(df) <- colnames(matrix)
  
  N <- length(unique(df$g)) # number of groups
  
  combinations <- data.frame(g1 = rep(unique(df$g), each = N), g2 = rep(unique(df$g), N), stringsAsFactors = FALSE)
  combinations <- combinations[combinations$g1 != combinations$g2, ]
  
  combinations$b1 <- df$b[match(combinations$g1, df$g)]
  combinations$b2 <- df$b[match(combinations$g2, df$g)]
  combinations <- combinations[combinations$b1!=combinations$b2,]
  
  idx <- c()
  for(i in 2:nrow(combinations)){
    if(!combinations$g2[i] %in% combinations$g1[1:(i-1)]) {
      idx <- c(idx, i)
    }
  }
  
  combinations <- combinations[c(1,idx),]
  rownames(combinations) <- 1:nrow(combinations)

  dist_p <- dist_t <- dist_coef <- matrix(0, nrow = N, ncol = N)
  colnames(dist_p) <- rownames(dist_p) <- unique(df$g)
  colnames(dist_t) <- rownames(dist_t) <- unique(df$g)
  colnames(dist_coef) <- rownames(dist_coef) <- unique(df$g)
  
  # create progress bar
  if(verbose == TRUE) {
    
    message("Generating distance matrix...")
    pb <- txtProgressBar(min = 0, max = nrow(combinations), style = 3)
    k <- 1
  }
  
  for(i in 1:nrow(combinations)){
    
    if(verbose == TRUE) {
      setTxtProgressBar(pb, k) # progress bar
      k <- k+1
    }
    
    df$tmp <- NA
    df$tmp <- "bg"
    df$tmp[df$g == combinations$g1[i]] <- "g1"
    df$tmp[df$g == combinations$g2[i]] <- "g2"
    df2 <- df[!is.na(df$tmp),]
    dge2 <- dge[,!is.na(df$tmp)]
    
    ## downsampling the bg
    n_bg <- sum(df2$tmp == "bg")
    
    if(bg.downsampling.factor > 1){
      set.seed(12345)
      random.idx <- sample(x = which(df2$tmp == "bg"), 
                           size = max(n_bg/bg.downsampling.factor , 50), replace = FALSE)
      select2 <- c(random.idx, which(df2$tmp %in% c("g1", "g2")) )
      select2 <- sort(select2)
    } else {
      select2 <- 1:n_bg
    }
    
    ## by group
    design <- model.matrix(~  0 + tmp + b + detrate, data = df2[select2,]) 
    groups <- paste0("tmp", unique(df2$tmp))
    groups <-  groups[groups != "tmpbg"]
    perm_groups <- data.frame(g1 = groups,
                              g2 = "tmpbg", stringsAsFactors = F)
    perm_groups <- perm_groups[perm_groups$g1 != perm_groups$g2,]
    perm_groups$pair <- paste0(perm_groups$g1, "-", perm_groups$g2)
    contrast_m <- makeContrasts(contrasts = perm_groups$pair,
                                levels = design)
    
    
    if(method == "voom"){
      v <- voom(dge2[,select2], design, plot = FALSE)
      fit <- lmFit(v, design)
      group_fit <- contrasts.fit(fit, contrast_m) 
      group_fit <- eBayes(group_fit)
      
    } else if (method == "trend") {
      logCPM <- cpm(dge2[,select2], log = TRUE, prior.count = 3)
      fit <- lmFit(logCPM, design)
      group_fit <- contrasts.fit(fit, contrast_m)
      group_fit <- eBayes(group_fit, trend = TRUE, robust = TRUE)
    }
 
    group_fit$p.value[,1] <- group_fit$p.value[,1] + 0.00000001
    group_fit$p.value[,2] <- group_fit$p.value[,2] + 0.00000001
    
    idx1 <- rownames(dist_coef) == combinations$g1[i]
    idx2 <- colnames(dist_coef) == combinations$g2[i]
    dist_coef[idx1, idx2] <- cor(coef(group_fit)[,1], coef(group_fit)[,2])
    dist_t   [idx1, idx2] <- cor(group_fit$t[,1], group_fit$t[,2])
    dist_p   [idx1, idx2] <- cor(-log10(group_fit$p.value)[,1]*sign(coef(group_fit)[,1]), 
                                 -log10(group_fit$p.value)[,2]*sign(coef(group_fit)[,2]))
    
  }
  
  if(verbose == TRUE) {
    close(pb) # close progress bar
  }
  
  return(list(dist_coef, dist_t, dist_p))
  
}


