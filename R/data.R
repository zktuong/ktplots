#' kidneyimmune
#'
#' A small set of demo data from Stewart et al. 2019 Science
#' @docType data
#' @usage data(kidneyimmune)
#' @format A Seurat object with the following slots filled
#' \describe{
#'   \item{assays}{
#'   \itemize{Currently only contains one assay ("RNA" - scRNA-seq expression data)
#'   \item{data - Normalized expression data}}
#' }
#'   \item{meta.data}{Cell level metadata}
#'   \item{active.assay}{Current default assay}
#'   \item{active.ident}{Current default idents}
#'   \item{version}{Seurat version used to create the object}
#' }
#' @source \url{https://www.kidneycellatlas.org/}
#' @examples
#' data(kidneyimmune)
"kidneyimmune"