#' @title Compute IDEr
#'
#' @description Calculate the similarity matrix based on the metrics of inter-group differential expression (IDEr)
#'
#' @param seu Seurat S4 object with the column of `initial_cluster` in its meta.data. Required.
#' @param group.by.var Character. Default: "initial_cluster".
#' @param method Character. It can be voom (default) or trend.
#' @param verbose Boolean. Print the message and progress bar. (Default: TRUE)
#' @param use.parallel Boolean. Use parallel. (Default: FALSE)
#' @param n.cores Numeric. Number of cores used for parallel computing. If no value is given (default), it will use the output of `detectCores(logical = FALSE)`.
#' @param downsampling.size Numeric. The number of cells representing each group. (Default: 40)
#' @param downsampling.include Boolean. Using `include = TRUE` to include the group smaller than required size. (Default: FALSE)
#' @param downsampling.replace Boolean. Using `replace = TRUE` if the group is smaller than required size and some cells will be repeatedly used. (Default: FALSE)
#' @param bg.downsampling.factor Numeric. The factor used to downsampling background cells. Using a number over 1 can decrease the computing time. The minimum number of background cells used is 50. (Default: 1)
#'
#' @details Details
#'
#' @return A list of four values: three similarity matrices and one list of indeces.
#'
#' @seealso \code{\link{plotNetwork}} \code{\link{finalClustering}}
#'
#' @export
#'
#' @import limma edgeR foreach utils doParallel
#' @importFrom parallel detectCores
#' @importFrom stats model.matrix cor
#'
getIDEr <- function(seu,
                    group.by.var = "initial_cluster",
                    batch.by.var = "Batch",
                    method = "voom",
                    verbose = TRUE,
                    use.parallel = FALSE,
                    n.cores = NULL,
                    downsampling.size = 40,
                    downsampling.include = TRUE,
                    downsampling.replace = TRUE,
                    bg.downsampling.factor = 1) {

  if(!group.by.var %in% colnames(seu@meta.data)) {
    warning("group.by.var is not in the colnames of seu@meta.data.")
    return(NULL)
  } 
  
  if(!batch.by.var %in% colnames(seu@meta.data)) {
    warning("batch.by.var is not in the colnames of seu@meta.data.")
    return(NULL)
  }
  
  tmp <- seu@meta.data[,colnames(seu@meta.data) == group.by.var]
  
  if(batch.by.var != "Batch") {
    seu$Batch <- seu@meta.data[,colnames(seu@meta.data) == batch.by.var]
  }
  
  ## merge seurat list
  metadata <- data.frame(
    label = tmp,
    batch = seu$Batch,
    # ground_truth = seu$Group, 
    stringsAsFactors = FALSE
  )

  select <- downsampling(
    metadata = metadata, n.size = downsampling.size,
    include = downsampling.include, replace = downsampling.replace
  )

  matrix <- as.matrix(seu@assays$RNA@counts[, select])
  colnames(matrix) <- paste0(colnames(matrix), 1:ncol(matrix)) # avoid duplication
  keep <- rowSums(matrix > 0.5) > 5
  dge <- edgeR::DGEList(counts = matrix[keep, , drop = F]) # make a edgeR object
  dge <- dge[!grepl("ERCC-", rownames(dge)), ] # remove ERCC
  dge <- dge[!grepl("MT-", rownames(dge)), ]

  df <- data.frame(
    g = metadata$label[select],
    b = metadata$batch[select], ## batch
    # ground_truth = metadata$ground_truth[select],
    stringsAsFactors = F
  ) ## label

  df$detrate <- scale(colMeans(matrix > 0))[, 1]
  rownames(df) <- colnames(matrix)

  N <- length(unique(df$g)) # number of groups

  # get the dataframe of combinations/pairs for comparison
  combinations <- data.frame(g1 = rep(unique(df$g), each = N), g2 = rep(unique(df$g), N), stringsAsFactors = FALSE)
  combinations <- combinations[combinations$g1 != combinations$g2, ]
  combinations$b1 <- df$b[match(combinations$g1, df$g)]
  combinations$b2 <- df$b[match(combinations$g2, df$g)]
  combinations <- combinations[combinations$b1 != combinations$b2, ]
  idx <- c()
  for (i in 2:nrow(combinations)) {
    if (!combinations$g2[i] %in% combinations$g1[1:(i - 1)]) {
      idx <- c(idx, i)
    }
  }
  combinations <- combinations[c(1, idx), ]
  rownames(combinations) <- 1:nrow(combinations)

  dist_p <- dist_t <- dist_coef <- matrix(0, nrow = N, ncol = N)
  colnames(dist_p) <- rownames(dist_p) <- unique(df$g)
  colnames(dist_t) <- rownames(dist_t) <- unique(df$g)
  colnames(dist_coef) <- rownames(dist_coef) <- unique(df$g)

  if (use.parallel == FALSE) {

    # create progress bar
    if (verbose == TRUE) {
      message("Generating distance matrix...")
      pb <- txtProgressBar(min = 0, max = nrow(combinations), style = 3)
      k <- 1
    }

    for (i in 1:nrow(combinations)) {
      if (verbose == TRUE) {
        setTxtProgressBar(pb, k) # progress bar
        k <- k + 1
      }

      df$tmp <- NA
      df$tmp <- "bg"
      df$tmp[df$g == combinations$g1[i]] <- "g1"
      df$tmp[df$g == combinations$g2[i]] <- "g2"
      df2 <- df[!is.na(df$tmp), ]
      dge2 <- dge[, !is.na(df$tmp)]

      ## downsampling the bg
      n_bg <- sum(df2$tmp == "bg")

      if (bg.downsampling.factor > 1) {
        set.seed(12345)
        random.idx <- sample(
          x = which(df2$tmp == "bg"),
          size = max(n_bg / bg.downsampling.factor, 50), replace = FALSE
        )
        select2 <- c(random.idx, which(df2$tmp %in% c("g1", "g2")))
        select2 <- sort(select2)
      } else {
        select2 <- 1:nrow(df2)
      }

      design <- model.matrix(~ 0 + tmp + b + detrate, data = df2[select2, ])
      groups <- paste0("tmp", unique(df2$tmp))
      groups <- groups[groups != "tmpbg"]
      perm_groups <- data.frame(
        g1 = groups,
        g2 = "tmpbg", stringsAsFactors = F
      )
      perm_groups <- perm_groups[perm_groups$g1 != perm_groups$g2, ]
      perm_groups$pair <- paste0(perm_groups$g1, "-", perm_groups$g2)
      contrast_m <- makeContrasts(
        contrasts = perm_groups$pair,
        levels = design
      )
      group_fit <- getGroupFit(dge2[, select2], design, contrast_m, method)

      idx1 <- rownames(dist_coef) == combinations$g1[i]
      idx2 <- colnames(dist_coef) == combinations$g2[i]
      dist_coef[idx1, idx2] <- cor(coef(group_fit)[, 1], coef(group_fit)[, 2])
      dist_t   [idx1, idx2] <- cor(group_fit$t[, 1], group_fit$t[, 2])
      dist_p   [idx1, idx2] <- cor(
        -log10(group_fit$p.value)[, 1] * sign(coef(group_fit)[, 1]),
        -log10(group_fit$p.value)[, 2] * sign(coef(group_fit)[, 2])
      )
    }
    if (verbose == TRUE) {
      close(pb) # close progress bar
    }
  } else if (use.parallel == TRUE) { # use parallel

    if (is.null(n.cores)) {
      n.cores <- detectCores(logical = FALSE)
    } else {
      n.cores <- min(n.cores, detectCores(logical = FALSE))
    }
    registerDoParallel(n.cores)
    n.iter <- nrow(combinations)

    i <- NULL
    j <- NULL

    df_dist <- foreach(i = combinations$g1, j = combinations$g2, df = rep(list(df), n.iter), dge = rep(list(dge), n.iter), bg.downsampling.factor = rep(bg.downsampling.factor, n.iter), method = rep(method, n.iter), .combine = "rbind") %dopar% {
      df$tmp <- NA
      df$tmp <- "bg"
      df$tmp[df$g == i] <- "g1"
      df$tmp[df$g == j] <- "g2"
      df2 <- df[!is.na(df$tmp), ]
      dge2 <- dge[, !is.na(df$tmp)]

      n_bg <- sum(df2$tmp == "bg")
      if (bg.downsampling.factor > 1) {
        set.seed(12345)
        random.idx <- sample(
          x = which(df2$tmp == "bg"),
          size = max(n_bg / bg.downsampling.factor, 50), replace = FALSE
        )
        select2 <- c(random.idx, which(df2$tmp %in% c("g1", "g2")))
        select2 <- sort(select2)
      } else {
        select2 <- 1:nrow(df2)
      }

      ## by group
      design <- model.matrix(~ 0 + tmp + b + detrate, data = df2[select2, ])
      contrast_m <- makeContrasts(
        contrasts = c("tmpg1-tmpbg", "tmpg2-tmpbg"),
        levels = design
      )
      group_fit <- getGroupFit(dge2[, select2], design, contrast_m, method)

      dist_coef <- cor(coef(group_fit)[, 1], coef(group_fit)[, 2])
      dist_t <- cor(group_fit$t[, 1], group_fit$t[, 2])
      dist_p <- cor(
        -log10(group_fit$p.value)[, 1] * sign(coef(group_fit)[, 1]),
        -log10(group_fit$p.value)[, 2] * sign(coef(group_fit)[, 2])
      )

      print(c(dist_coef, dist_t, dist_p))
    }


    for (i in 1:nrow(combinations)) {
      idx1 <- rownames(dist_coef) == combinations$g1[i]
      idx2 <- colnames(dist_coef) == combinations$g2[i]
      dist_coef[idx1, idx2] <- df_dist[i, 1]
      dist_t   [idx1, idx2] <- df_dist[i, 2]
      dist_p   [idx1, idx2] <- df_dist[i, 3]
    }
  }

  return(list(dist_coef, dist_t, dist_p, select, combinations))
}


#' @title Calculate Fit
#'
#' @param dge dge
#' @param design design
#' @param contrast_m contrast_m
#' @param method method
#'
#' @import limma edgeR
#' @export
#'
#' @seealso \code{\link{getIDEr}}
#'
getGroupFit <- function(dge, design, contrast_m, method) {
  if (method == "voom") {
    v <- voom(dge, design, plot = FALSE)
    fit <- lmFit(v, design)
    group_fit <- contrasts.fit(fit, contrast_m)
    group_fit <- eBayes(group_fit)
  } else if (method == "trend") {
    logCPM <- cpm(dge, log = TRUE, prior.count = 3)
    fit <- lmFit(logCPM, design)
    group_fit <- contrasts.fit(fit, contrast_m)
    group_fit <- eBayes(group_fit, trend = TRUE, robust = TRUE)
  }

  group_fit$p.value[, 1] <- group_fit$p.value[, 1] + 0.00000001
  group_fit$p.value[, 2] <- group_fit$p.value[, 2] + 0.00000001

  return(group_fit)
}
