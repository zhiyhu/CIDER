#' @title Initial clustering for evaluating integration
#'
#' @description This function applies HDBSCAN, a density-based
#' clustering method, on the corrected dimension reduction.
#'
#' @param seu a Seurat object containing integrated or batch corrected
#'  PCA.
#' @param reduction Character. Name of the dimension reduction after
#' integration or batch correction. (Default: PCA)
#' @param dims Numeric vector. Dimensions used for initial clustering.
#' (Default: 1:15)
#' @param minPts Interger. Minimum size of clusters. Will be passed
#' to the `hdbscan` function. (Default: 25)
#'
#' @return A Seurat object having two additional columns in its
#' meta.data: dbscan_cluster and initial_cluster.
#'
#' @seealso Usage of this function should be followed by
#' \code{\link{getIDEr}} and \code{\link{estimateProb}}.
#' @export
#'
#' @import Seurat
#' @importFrom dbscan hdbscan
hdbscan.seurat <- function(seu, reduction = "pca",
                            dims = seq_len(15), minPts = 25){
  if(!reduction %in% Reductions(seu)) stop("pca does not exist.")
  seu <- RunTSNE(seu, reduction = reduction, dims = dims)
  res <- hdbscan(Reductions(seu, "tsne")@cell.embeddings[,seq_len(2)],
                  minPts = minPts)
  seu$dbscan_cluster <- factor(as.character(res$cluster))
  seu$initial_cluster <- factor(as.character(paste0(seu$dbscan_cluster,
                                                    "_", seu$Batch)))
  return(seu)
}


#' @title Estimate the empirical probability of whether two set of cells
#' from distinct batches belong to the same population
#' @param seu A Seurat object
#' @param ider The output list of function `getIDEr`.
#' @param n_size Number of cells per group used to compute the similarity. Default: 40
#' @param n.perm Numeric. Time of permutations.
#' @param verbose Boolean. Print out progress or not. (Default: FALSEW)
#' @return A Seurat object with IDER-based similarity and empirical
#' probability of rejection
#' @import limma edgeR foreach utils doParallel
#' @importFrom kernlab specc
#' @seealso Usage of this function should be after \code{\link{hdbscan.seurat}}
#' and \code{\link{getIDEr}}
#' @export
estimateProb <- function(seu, ider, n_size = 40,
                          #seeds = c(12345, 89465, 10385, 10385, 17396),
                          n.perm = 5, verbose = FALSE){

  dist_coef <- ider[[1]]
  dist_coef[upper.tri(dist_coef)] <- 0
  # Select positive control
  pos_control <- c(rownames(dist_coef)[which.max(apply(dist_coef, 1, max))],
                   colnames(dist_coef)[which.max(apply(dist_coef, 2, max))])

  idx <- seu$initial_cluster %in% pos_control

  combinations_all <- c()
  bg_dist_coef_list <- list() # background distance distribution
  seu_selected <- seu[,idx]
  pca <- seu@reductions$pca@cell.embeddings[idx, seq_len(15)]
  # use first 15 PCs for spectral clustering

  for(itor in seq_len(n.perm)) {
    # force the positive control group into first groups
    # set.seed(seeds[itor])
    res <- specc(pca, centers = 5) # spectral clustering

    ## Calculate background distribution
    seu_selected$forced_cluster <- res@.Data
    seu$forced_cluster <-
      seu_selected$forced_cluster[match(colnames(seu),
                                        colnames(seu_selected))]
    seu$forced_cluster[is.na(seu$forced_cluster)] <- "bg"
    metadata <- data.frame(label = paste0(seu$forced_cluster, "_", seu$Batch),
                            batch = seu$Batch,
                            stringsAsFactors = FALSE)

    # downsampling
    select <- downsampling(metadata = metadata, n.size = n_size,
                           include = TRUE, replace = TRUE, lower.cutoff = 15)

    # IDER
    matrix <- as.matrix(seu@assays$RNA@counts[,select])
    keep <- rowSums(matrix > 0.5) > 5
    dge <- edgeR::DGEList(counts = matrix[keep, , drop = FALSE])
    # make a edgeR object
    dge <- dge[!grepl("ERCC-", rownames(dge)),] # remove ERCC
    dge <- dge[!grepl("MT-", rownames(dge)),]

    df <- data.frame(g = metadata$label[select],
                     b = metadata$batch[select], ## batch
                     stringsAsFactors = FALSE) ## label
    df$detrate <- scale(colMeans(matrix > 0))[, 1]
    colnames(matrix) <- paste0(colnames(matrix), "_", seq_len(ncol(matrix)))
    rownames(df) <- colnames(matrix)

    GROUPS <- unique(df$g)
    N <- length(GROUPS)

    combinations <- data.frame(g1 = rep(unique(df$g), each = N),
                               g2 = rep(unique(df$g), N),
                               stringsAsFactors = FALSE)
    combinations <- combinations[combinations$g1 != combinations$g2, ]
    combinations$b1 <- df$b[match(combinations$g1, df$g)]
    combinations$b2 <- df$b[match(combinations$g2, df$g)]
    combinations <- combinations[combinations$b1 != combinations$b2, ]
    idx <- c()
    for(i in 2:nrow(combinations)){
      if(!combinations$g2[i] %in% combinations$g1[seq_len(i-1)]) {
        idx <- c(idx, i)
      }
    }

    combinations <- combinations[c(1,idx),]
    rownames(combinations) <- seq_len(nrow(combinations))
    combinations <-
      combinations[!combinations$g1 %in% c("bg_Batch1", "bg_Batch2")  &
                     !combinations$g2 %in% c("bg_Batch1", "bg_Batch2"),]
    combinations$similarity <- NA
    combinations$iteration <- itor

    bg_dist_coef <- matrix(0, nrow = N, ncol = N)
    colnames(bg_dist_coef) <- rownames(bg_dist_coef) <- GROUPS

    # create progress bar
    if (verbose == TRUE) {
      message("Generating distance matrix...")
      pb <- txtProgressBar(min = 0, max = nrow(combinations), style = 3)
      k <- 1
    }
    logCPM_all <- cpm(dge, log = TRUE, prior.count = 3)
    for (i in seq_len(nrow(combinations))){
      if (verbose == TRUE) {
        setTxtProgressBar(pb, k) # progress bar
        k <- k+1
      }
      df$tmp <- NA
      df$tmp[df$g %in% c("bg_Batch1", "bg_Batch2")] <- "bg"
      df$tmp[df$g == combinations$g1[i]] <- "g1"
      df$tmp[df$g == combinations$g2[i]] <- "g2"

      idx <- which(!is.na(df$tmp))
      design <- model.matrix(~  0 + tmp + b + detrate, data = df[idx, ])
      contrast_m <- makeContrasts(contrasts = c("tmpg1-tmpbg", "tmpg2-tmpbg"),
                                  levels = design)
      logCPM <- logCPM_all[,idx]
      fit <- lmFit(logCPM, design)
      group_fit <- contrasts.fit(fit, contrast_m)

      idx1 <- rownames(bg_dist_coef) == combinations$g1[i]
      idx2 <- colnames(bg_dist_coef) == combinations$g2[i]
      combinations$similarity[i] <- bg_dist_coef[idx1, idx2] <-
        cor(coef(group_fit)[,1], coef(group_fit)[,2])
    }

    combinations_all <- rbind(combinations_all, combinations)
    bg_dist_coef_list[[itor]] <- bg_dist_coef

    if(verbose == TRUE) {
      close(pb) # close progress bar
    }
  }

  idx <- getSharedGroups(seu, ider[[1]])
  shared_g <- idx[[1]]
  idx1 <- idx[[2]]
  idx2 <- idx[[3]]

  p_mat <- getProbability(ider[[1]][idx1, idx2], combinations_all$similarity)

  # assign similiary
  scores <- diag(ider[[1]][idx1, idx2])
  names(scores) <- shared_g
  seu$similarity <- scores[match(seu$dbscan_cluster, names(scores))]

  # assign p values
  scores <- diag(p_mat[idx1, idx2])
  names(scores) <- shared_g
  seu$pvalue <- scores[match(seu$dbscan_cluster, names(scores))]

  return(seu)
}

getProbability <- function(x, bg_similarity) { # calculate the probability
  try({
    res <- matrix(NA, ncol = ncol(x), nrow = nrow(x))
    rownames(res) <- rownames(x)
    colnames(res) <- colnames(x)
    len <- length(bg_similarity)
    for(i in seq_len(nrow(x))){
      for(j in seq_len(ncol(x))){
        res[i,j] <- sum(bg_similarity > x[i,j]) / len
      }
    }
    return(res)
  })
}
