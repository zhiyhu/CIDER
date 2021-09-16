#' Count matrix of pancreatic scRNA-Seq data.
#'
#' A count matrix (dgCMatrix) containing pancreatic scRNA-Seq data from both human and mouse.
#'
#' @format A sparse matrix with xxx rows and xxx variables:
#' \describe{
#'   \item{Rows}{Gene names.}
#'   \item{Col2}{Cell names.}
#' }
#' @source \url{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE84133}
#' @examples 
#' data("pancreas_meta")
"pancreas_counts"

#' Metadata of pancreatic scRNA-Seq data.
#'
#' A data.frame containing cell annotations (meta data).
#' 
#' @format A data frame with xxx rows and xxx variables:
#' \describe{
#'   \item{Col1}{Batch}
#'   \item{Col2}{Group}
#'   \item{Col3}{Sample}
#' }
#' @source \url{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE84133}
#' @examples 
#' data("pancreas_meta")
"pancreas_meta"

#' Seurat object of dendritic scRNA-Seq data.
#'
#' A Seurat object of dendritic scRNA-Seq data.
#' 
#' @format A Seurat object Containing xxxxx.
#' @source \url{https://www.ncbi.nlm.nih.gov/geo}
#' @examples 
#' data("dendritic")
"dendritic"