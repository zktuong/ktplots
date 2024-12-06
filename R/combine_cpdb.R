#' combine multiple results from cellphonedb.
#'
#' @param ... multiple cellphonedb outputs.
#' @return combine results from multiple cellphonedb runs.
#' @import dplyr
#' @import purrr
#' @export
combine_cpdb <- function(...) {
    output <- list(...)
    anames <- c(
        "id_cp_interaction", "interacting_pair", "partner_a", "partner_b",
        "gene_a", "gene_b", "secreted", "receptor_a", "receptor_b", "annotation_strategy",
        "is_integrin"
    )
    bnames <- c(
        "gene_name", "uniprot", "is_complex", "protein_name", "complex_name",
        "id_cp_interaction"
    )
    if (all(colnames(output[[1]])[1:11] == anames)) {
        out <- output %>%
            reduce(full_join, by = anames)
    } else if (all(colnames(output[[1]])[1:6] == bnames)) {
        out <- output %>%
            reduce(full_join, by = bnames)
    }
    return(out)
}
