#' kidneyimmune
#'
#' kidneyimmune - A small set of demo data from Stewart et al. 2019 Science. See \url{https://www.kidneycellatlas.org/}
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

#' means
#'
#' means - Dataframe of CellPhoneDB output means.txt file
#' @rdname kidneyimmune
#' @docType data
#' @usage data(cpdb_output)
#' @format Data frames containing outputs after CellPhoneDB analysis
#' \describe{
#'   \item{means}{data frame containing mean expression values for each interacting pair}
#'   \item{pvals}{data frame containing p values for each interacting pair}
#' }
#' @examples
#' data(cpdb_output)
"means"

#' pvals
#'
#' pvals - Dataframe of CellPhoneDB output pvalues.txt file
#' @rdname kidneyimmune
#' @docType data
"pvals"


#' decon
#'
#' decon - Dataframe of CellPhoneDB output deconvoluted.txt file
#' @rdname kidneyimmune
#' @docType data
"decon"


#' sig_means
#'
#' sig_means - Dataframe of CellPhoneDB output significant_means.txt file
#' @rdname kidneyimmune
#' @docType data
"sig_means"


#' interaction_annotation
#'
#' interaction_annotation - Example dataframe to use for interaction_grouping option in plot_cpdb2
#' @rdname kidneyimmune
#' @docType data
"interaction_annotation"


