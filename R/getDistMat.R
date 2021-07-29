#' @title Calculate Distance Matrix
#'
#' @description discription
#'
#' @author Zhiyuan Hu
#'
#' @param seu_list seu list
#' @param verbose print the message and progress bar (default: TRUE)
#' @param tmp.initial.clusters tmp.initial.clusters
#' @param method voom or trend
#' @param additional.variate additional variate to include into the regression model
#' @param downsampling.size downsampling.size
#' @param downsampling.include downsampling.include
#' @param downsampling.replace downsampling.replace
#'
#' @return a distance matrix
#'
#' @export
#'
#' @import Seurat utils limma
#' @importFrom edgeR cpm
#'
getDistMat <- function(seu_list,
                       verbose = TRUE,
                       tmp.initial.clusters = "seurat_clusters",
                       method = "trend",
                       additional.variate = NULL,
                       downsampling.size = 35,
                       downsampling.include = TRUE,
                       downsampling.replace = TRUE) {
  dist_coef <- list()

  if (verbose == TRUE) { # create progress bar
    pb <- txtProgressBar(min = 0, max = length(seu_list), style = 3)
    k <- 1
  }

  for (seu_itor in 1:length(seu_list)) {
    df_info <- data.frame(
      label = seu_list[[seu_itor]]$seurat_clusters,
      batch = seu_list[[seu_itor]]$Batch
      # donor = seu_list[[seu_itor]]$Tissue
    )

    idx <- downsampling(
      metadata = df_info, n.size = downsampling.size,
      include = downsampling.include, replace = downsampling.replace
    )
    idx <- sort(idx)

    to_add <- idx[duplicated(idx)]
    idx <- idx[!duplicated(idx)]
    matrix <- as.matrix(seu_list[[seu_itor]]@assays$RNA@counts[, idx])

    if (length(to_add) > 0) {
      matrix2 <- data.frame(seu_list[[seu_itor]]@assays$RNA@counts[, to_add])
      colnames(matrix2) <- paste0(colnames(matrix2), 1:ncol(matrix2))
      matrix2 <- as.matrix(matrix2)
      matrix <- cbind(matrix, matrix2)
      rm(matrix2)
    }

    if (length(unique(df_info$label[idx])) > 2) {
      dist_coef[[seu_itor]] <- calculateDistMatOneModel(
        matrix = matrix, metadata = df_info[c(idx, to_add), ],
        # matrix = matrix, metadata = df_info[idx, ],
        verbose = verbose, method = method,
        additional.variate = additional.variate
      )
    }
    
    if (verbose == TRUE) {
      setTxtProgressBar(pb, k) # progress bar
      k <- k + 1
    }
  }
  if (verbose == TRUE) {
    close(pb) # close progress bar
  }

  return(dist_coef)
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
#' @import limma edgeR
#' @importFrom stats model.matrix cor coef
#'
calculateDistMatOneModel <- function(matrix, metadata,
                                     verbose = TRUE,
                                     method = "voom",
                                     additional.variate = NULL) 
  {
  keep <- rowSums(matrix > 0.5) > 5 
  dge <- edgeR::DGEList(counts = matrix[keep,,drop=F]) # make a edgeR object
  dge <- dge[!grepl("ERCC-", rownames(dge)),] # remove ERCC
  dge <- dge[!grepl("MT-", rownames(dge)),]
  dge <- dge[!grepl("mt-", rownames(dge)),]

  df <- data.frame(g = paste(metadata$label, metadata$batch, sep = "_"),
                   b = metadata$batch, ## batch
                   c = metadata$label, ## label
                   stringsAsFactors = F) 
  df$detrate <- scale(colMeans(matrix > 0))[,1] # gene detection rate
  rownames(df) <- colnames(matrix)
  
  N <- length(unique(df$g)) # number of initial groups
  
  combinations <- data.frame(g1 = rep(unique(df$g), each = N), g2 = rep(unique(df$g), N), stringsAsFactors = FALSE)
  combinations <- combinations[combinations$g1 != combinations$g2, ]
  # combinations <- sort(combinations)
  
  idx <- c()
  for(i in 2:nrow(combinations)){
    if(!combinations$g2[i] %in% combinations$g1[1:(i-1)]) {
      idx <- c(idx, i)
    }
  }
  
  combinations <- combinations[c(1,idx),]
  rownames(combinations) <- 1:nrow(combinations)
  
  dist_coef <- matrix(0, nrow = N, ncol = N)
  colnames(dist_coef) <- rownames(dist_coef) <- sort(unique(df$g))
  
  if("donor" %in% colnames(metadata)){
    df$subb <- metadata$donor
    design <- model.matrix(~  0 + g + subb + detrate, data = df)
  } else {
    design <- model.matrix(~  0 + g + detrate, data = df)
  }
  
  groups <- sort(unique(paste0("g", df$g)))
  n_groups <- length(groups) # number of groups
  
  df_contrasts <- data.frame(target_group = groups, contrast = NA) # prepare contrast matrix
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
    # group_fit <- eBayes(group_fit)
  } else if (method == "trend") {
    logCPM <- edgeR::cpm(dge, log=TRUE, prior.count=3)
    fit <- lmFit(logCPM, design)
    group_fit <- contrasts.fit(fit, contrast_m)
    # group_fit <- eBayes(fit, trend=TRUE, robust = TRUE)
  }
  # pairwise comparison
  for(i in 1:nrow(combinations)){
    idx1 <- rownames(dist_coef) == combinations$g1[i]
    idx2 <- colnames(dist_coef) == combinations$g2[i]
    pos1 <- df_contrasts$target_group == paste0("g", combinations$g1[i])
    pos2 <- df_contrasts$target_group == paste0("g", combinations$g2[i])
    
    dist_coef[idx1, idx2] <- cor(coef(group_fit)[, pos1], coef(group_fit)[, pos2])
  }
  return(dist_coef)
}

