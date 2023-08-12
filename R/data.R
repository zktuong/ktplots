#' kidneyimmune
#'
#' kidneyimmune - A small set of demo data from Stewart et al. 2019 Science. See \url{https://www.kidneycellatlas.org/}
#' @docType data
#' @usage data(kidneyimmune)
#' @format A SingleCellExperiment object with the following slots filled
#' \describe{
#'   \item{assays}{
#'   \itemize{Currently only contains "counts" and "logcounts"
#'   \item{counts - Raw expression data}
#'   \item{logcounts - Normalized expression data}}
#' }
#'   \item{colData}{Cell level metadata}
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

#' means2
#'
#' means2 - Dataframe of CellPhoneDB output means.txt file
#' @rdname kidneyimmune
#' @docType data
#' @usage data(cpdb_output)
#' @format Data frames containing outputs after CellPhoneDB analysis
#' \describe{
#'   \item{means2}{data frame containing mean expression values for each interacting pair}
#'   \item{pvals2}{data frame containing p values for each interacting pair}
#' }
#' @examples
#' data(cpdb_output)
"means2"

#' pvals2
#'
#' pvals2 - Dataframe of CellPhoneDB output pvalues.txt file
#' @rdname kidneyimmune
#' @docType data
"pvals2"


#' decon2
#'
#' decon2 - Dataframe of CellPhoneDB output deconvoluted.txt file
#' @rdname kidneyimmune
#' @docType data
"decon2"


#' sig_means2
#'
#' sig_means2 - Dataframe of CellPhoneDB output significant_means.txt file
#' @rdname kidneyimmune
#' @docType data
"sig_means2"


#' covid_cpdb_meta
#'
#' covid_cpdb_meta - Example dataframe to use for cpdb_meta option in compare_cpdb
#' @rdname kidneyimmune
#' @docType data
"covid_cpdb_meta"


#' covid_sample_metadata
#'
#' covid_sample_metadata - Example dataframe to use for sample_metadata option in compare_cpdb
#' @rdname kidneyimmune
#' @docType data
"covid_sample_metadata"

#' cpdb_output_v5
#'
#' cpdb_output_v5 - Dataframe of CellPhoneDB output means.txt file
#' @rdname kidneyimmune
#' @docType data
#' @usage data(cpdb_output_v5)
#' @format data after CellPhoneDB v5 analysis
#' @examples
#' data(cpdb_output_v5)
"means_v5"


#' sce_v5
#'
#' sce_v5 - A small dummy singlecelldata for cellphonedb v5
#' @rdname kidneyimmune
#' @docType data
"sce_v5"


#' relevant_interactions_v5
#'
#' relevant_interactions_v5 - Dataframe of CellPhoneDB output relevant_interactions.txt file
#' @rdname kidneyimmune
#' @docType data
"relevant_interactions_v5"


#' decon_v5
#'
#' decon_v5 - Dataframe of CellPhoneDB output deconvoluted.txt file
#' @rdname kidneyimmune
#' @docType data
"decon_v5"


#' cellsign_v5
#'
#' cellsign_v5 - Dataframe of CellPhoneDB output CellSign.txt file
#' @rdname kidneyimmune
#' @docType data
"cellsign_v5"

#' interaction_scores_v5
#'
#' interaction_scores_v5 - Dataframe of CellPhoneDB output interaction_scores.txt file
#' @rdname kidneyimmune
#' @docType data
"interaction_scores_v5"
