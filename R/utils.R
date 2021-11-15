#' @title Downsampling cells
#'
#' @description Downsampling cells from each group for IDER-based similarity
#' calculation.
#'
#' @param metadata Data frame. It includes at least 2 columns, label and batch.
#' Each row corresponds to one cell. Required.
#' @param n.size Numeric. The number of cells used in each group. (Default: 35)
#' @param seed Numeric. Seed used to sample. (Default: 12345)
#' @param include Boolean. Using `include = TRUE` to include the group smaller
#' than required size. (Default: FALSE)
#' @param replace Boolean. Using `replace = TRUE` if the group is smaller than
#' required size and some cells will be repeatedly used. (Default: FALSE)
#' @param lower.cutoff Numeric. The minimum size of groups to keep.
#' (Default: 3)
#'
#' @return A numeric list of which cells will be kept for downstream
#' computation.
#'
#' @export
#'
downsampling <- function(metadata, n.size = 35, seed = 12345, include = FALSE,
                         replace = FALSE, lower.cutoff = 3) {
  cluster <- unique(metadata$label)
  tech <- unique(metadata$batch)
  select <- c()
  for (i in cluster) {
    for (j in tech) {
      idx <- which(metadata$label %in% i & metadata$batch %in% j)
      if (length(idx) > n.size) {
        # set.seed(seed)
        select <- c(select, sample(idx, size = n.size, replace = FALSE))
      } else if (length(idx) == n.size) {
        select <- idx
      } else if (length(idx) < n.size & length(idx) >= lower.cutoff) {
        if (include & !replace) {
          select <- c(select, idx)
        } else if (include & replace) {
          # set.seed(seed)
          select <- c(select, sample(idx, size = n.size, replace = TRUE))
        }
      }
    }
  }
  return(select)
}

dist2similarity <- function(dist){
  sml <- dist[[1]] + t(dist[[1]])
  sml <- 1 - sml
  diag(sml) <- 1
}

getSharedGroups <- function(seu, dist){

  batches <- unique(seu$Batch)
  groups <- colnames(dist)
  idx1 <- colnames(dist)[grep(paste0("_", batches[1],"$"), colnames(dist))]
  # end with _Batch1
  idx2 <- colnames(dist)[grep(paste0("_", batches[2],"$"), colnames(dist))]
  # end with _Batch2

  g1 <- gsub(paste0("_", batches[1]), "",idx1)
  g2 <- gsub(paste0("_", batches[2]), "",idx2)
  shared_g <- intersect(g1,g2)

  idx1_shared <- paste0(shared_g, paste0("_", batches[1]))
  idx2_shared <- paste0(shared_g, paste0("_", batches[2]))

  return(list(shared_g, idx1_shared, idx2_shared))
}


#' @title Measure similarity between two vectors
#' @description Measure similarity between two vectors
#' @param x1 x1
#' @param x2 x2
#' @param method method
#' @return similarity matrix
#' @importFrom stats cor
measureSimilarity <- function(x1, x2, method = "pearson"){
  if (!is.null(x1) & !is.null(x2)) {
    if(length(x1) != length(x2)) {
      warning("x1 or x2 don't have the same length for similarity measures")
      return(NA)
    }
    if (method %in% c("pearson", "spearman","kendall")) {
      return(cor(x1, x2, method = method))
    } else if (method == "cosine") {
      x <- matrix(cbind(x1, x2))
      return(1-cosineSimilarityR(x)[1,2])
    }
  } else {
    warning("x1 or x2 are not valid for similarity measures")
    return(NA)
  }
}

#' cosine similarity in R
#' @param x a matrix
#' @return a similarity matrix among all rows of the input matrix
cosineSimilarityR <- function(x) {
  y <- t(x) %*% x
  res <- y / (sqrt(diag(y)) %*% t(sqrt(diag(y))))
  return(res)
}


#' @references Ritchie ME, Phipson B, Wu D, Hu Y, Law CW, Shi W, Smyth GK
#' (2015). “limma powers differential expression analyses for RNA-sequencing
#' and microarray studies.” Nucleic Acids Research, 43(7), e47.
#' doi: 10.1093/nar/gkv007.
.zeroDominantMatrixMult <- function(A,B)
{
  # Computes A %*% B, except that a zero in B will always produce
  # zero even when multiplied by an NA in A, instead of NA as usually
  # produced by R arithmetic.
  # A and B are numeric matrices and B does not contain NAs.
  # In the limma usage, A usually has far more rows than columns
  # and B is relatively small.
  # AUTHOR: Gordon Smyth
  # Created 16 Feb 2018. Modified 2 Feb 2020.
  # THIS IS AN UNEXPORTED FUNCTION FROM R PACKAGE LIMMA
  # https://bioconductor.org/packages/release/bioc/html/limma.html

  # Proportion of zeros in B
  Z <- (B==0)
  MeanZ <- mean(Z)

  # Decide whether to run guarded or ordinary matrix multiplication
  # MeanZ can only be NA if B has 0 elements
  if(!is.na(MeanZ) && (MeanZ > 0)) {
    if(MeanZ >= 0.8)
      # Full algorithm is quick if there are lots of zeros
      Guard <- TRUE
    else {
      RowBHasZero <- (rowSums(Z) > 0)
      if(mean(RowBHasZero) > 0.4) {
        # If the matrix is big, it's much quicker to check the whole matrix
        # than to subset it
        Guard <- anyNA(A)
      } else {
        Guard <- anyNA(A[,RowBHasZero])
      }
    }
  } else {
    Guard <- FALSE
  }

  if(Guard) {
    dn <- list()
    dn[[1]] <- rownames(A)
    dn[[2]] <- colnames(B)
    D <- matrix(0,nrow(A),ncol(B),dimnames=dn)
    for (j in seq_len(ncol(B))) {
      z <- B[,j]==0
      if(any(z))
        D[,j] <- A[,!z,drop=FALSE] %*% B[!z,j,drop=FALSE]
      else
        D[,j] <- A %*% B[,j]
    }
    return(D)
  } else {
    return(A %*% B)
  }
}
