#' Plotting CellPhoneDB results as a chord diagram
#'
#' @param scdata single-cell data. Must be a SingleCellExperiment object.
#' @param cell_type1 cell type 1
#' @param cell_type2 cell type 2
#' @param celltype_key vector holding the celltype_key for each cell or column name of scdata's metadata. MUST match cpdb's columns
#' @param means object holding means.txt from cpdb output
#' @param pvals object holding pvals.txt from cpdb output
#' @param deconvoluted object holding deconvoluted.txt from cpdb output
#' @param keep_significant_only logical. Default is FALSE. Switch to TRUE if you only want to plot the significant hits from cpdb.
#' @param splitby_key column name in the metadata/coldata table to split the spots by. Can only take columns with binary options. If specified, name to split by MUST be specified in the meta file provided to cpdb prior to analysis.
#' @param standard_scale logical. scale the expression to range from 0 to 1. Default is TRUE
#' @param gene_symbol_mapping default = NULL.column name for rowData in sce holding the actual gene symbols if row names aren't gene symbols
#' @param frac default = 0.1. Minimum fraction of celtypes expressing a gene in order to keep the interaction. Gene must be expressesd >= `frac` in either of the pair of celltypes in order to keep.
#' @param remove_self default = TRUE. Remove self-self arcs.
#' @param desiredInteractions default = NULL. Specific list of celltype comparisons e.g. list(c('CD4_Tcm', 'cDC1'), c('CD4_Tcm', 'cDC2')). Also accepts a dataframe where first column is celltype 1 and 2nd column is celltype 2.
#' @param degs_analysis if is CellPhoneDB degs_analysis mode.
#' @param directional Whether links have directions. 1 means the direction is from the first column in df to the second column, -1 is the reverse, 0 is no direction, and 2 for two directional.
#' @param alpha transparency for links
#' @param edge_colors vector of colors for links
#' @param grid_colors vector of colors for grids
#' @param show_legend whether or not to show the legend
#' @param legend.pos.x x position of legend
#' @param legend.pos.y y position of legend
#' @param ... passes arguments plot_cpdb
#' @return Plotting CellPhoneDB results as a CellChat inspired chord diagram
#' @examples
#' \donttest{
#'
#' }
#' @include utils.R
#' @importFrom circlize circos.clear chordDiagram
#' @importFrom grDevices recordPlot
#' @export

plot_cpdb3 <- function(
    scdata, cell_type1, cell_type2, celltype_key, means, pvals,
    deconvoluted, keep_significant_only = TRUE, splitby_key = NULL, standard_scale = TRUE,
    gene_symbol_mapping = NULL, frac = 0.1, remove_self = TRUE, desiredInteractions = NULL,
    degs_analysis = FALSE, directional = 1, alpha = 0.5, edge_colors = NULL, grid_colors = NULL,
    show_legend = TRUE, legend.pos.x = 20, legend.pos.y = 20, ...) {
    if (class(scdata) == "Seurat") {
        stop("Sorry not supported. Please use a SingleCellExperiment object.")
    }
    lr_interactions <- plot_cpdb(
        scdata = scdata, cell_type1 = cell_type1, cell_type2 = cell_type2,
        celltype_key = celltype_key, splitby_key = splitby_key, means = means, pvals = pvals,
        keep_significant_only = keep_significant_only, standard_scale = standard_scale,
        return_table = TRUE, degs_analysis = degs_analysis, ...
    )
    if (class(lr_interactions) == "logical") { # check if NA
        message("No interactions could be plotted.")
        return(NA)
    }
    requireNamespace("SummarizedExperiment")
    requireNamespace("SingleCellExperiment")
    means_col <- grep("scaled_means|means", names(lr_interactions), value = TRUE)[1]
    if (is.null(splitby_key)) {
        if (any(lr_interactions[, means_col] > 0)) {
            if (any(is.na(lr_interactions[, means_col]))) {
                lr_interactions <- lr_interactions[lr_interactions[, means_col] >
                    0 & !is.na(lr_interactions[, means_col]), ]
            } else {
                lr_interactions <- lr_interactions[lr_interactions[, means_col] >
                    0, ]
            }
        }
    }
    subset_clusters <- unique(unlist(lapply(
        as.character(lr_interactions$group),
        strsplit, DEFAULT_SEP
    )))
    sce_subset <- scdata[, SummarizedExperiment::colData(scdata)[, celltype_key] %in%
        subset_clusters]
    interactions <- means[, c(
        "id_cp_interaction", "interacting_pair", "gene_a",
        "gene_b", "partner_a", "partner_b", "receptor_a", "receptor_b"
    )]
    interactions$use_interaction_name <- paste0(
        interactions$id_cp_interaction, SPECIAL_SEP,
        interactions$interacting_pair
    )
    interactions$converted <- gsub("_", "-", interactions$use_interaction_name)
    interactions_subset <- interactions[interactions$converted %in% lr_interactions$Var1, ]
    tm0 <- do.call(c, lapply(
        as.list(interactions_subset$use_interaction_name), strsplit,
        "_"
    ))
    if (any(lapply(tm0, length) > 2)) {
        complex_id <- which(lapply(tm0, length) > 2)
        interactions_subset_ <- interactions_subset[complex_id, ]
        simple_1 <- interactions_subset_$interacting_pair[grep("complex:", interactions_subset_$partner_b)]
        partner_1 <- gsub("complex:", "", interactions_subset_$partner_b[grep(
            "complex:",
            interactions_subset_$partner_b
        )])
        partner_2 <- gsub("complex:", "", interactions_subset_$partner_a[grep(
            "complex:",
            interactions_subset_$partner_a
        )])
        simple_2 <- interactions_subset_$interacting_pair[grep("complex:", interactions_subset_$partner_a)]
        for (i in seq_along(simple_1)) {
            simple_1[i] <- gsub(paste0(partner_1[i], "_|_", partner_1[i]), "", simple_1[i])
        }
        for (i in seq_along(simple_2)) {
            simple_2[i] <- gsub(paste0(partner_2[i], "_|_", partner_2[i]), "", simple_2[i])
        }
        tmpdf <- rbind(cbind(simple_1, partner_1), cbind(partner_2, simple_2))
        tmplist <- split(tmpdf, seq(nrow(tmpdf)))
        for (i in seq_along(complex_id)) {
            tm0[[complex_id[i]]] <- tmplist[[i]]
        }
        tm0 <- data.frame(t(matrix(unlist(tm0), 2, length(unlist(tm0)) / 2)))
        colnames(tm0) <- c("id_a", "id_b")
        interactions_subset <- cbind(interactions_subset, tm0)
        dictionary <- interactions_subset[, c(
            "id_cp_interaction", "gene_a", "gene_b",
            "partner_a", "partner_b", "id_a", "id_b", "receptor_a", "receptor_b"
        )]
    } else {
        tm0 <- data.frame(t(matrix(unlist(tm0), 2, length(unlist(tm0)) / 2)))
        colnames(tm0) <- c("id_a", "id_b")
        tm0$id_a <- gsub(paste0(".*", SPECIAL_SEP), "", tm0$id_a)
        interactions_subset <- cbind(interactions_subset, tm0)
        dictionary <- interactions_subset[, c(
            "id_cp_interaction", "gene_a", "gene_b",
            "partner_a", "partner_b", "id_a", "id_b", "receptor_a", "receptor_b"
        )]
    }
    # extract all the possible genes.
    geneid <- unique(c(interactions_subset$id_a, interactions_subset$id_b))
    if (all(!geneid %in% row.names(sce_subset))) {
        geneid <- unique(c(interactions_subset$gene_a, interactions_subset$gene_b))
    }
    sce_subset_tmp <- sce_subset[row.names(sce_subset) %in% geneid, ]
    if (dim(sce_subset_tmp)[1] == 0) {
        stop("Gene ids not found. Are you sure your single-cell object contains human gene symbols?")
    }
    # split to list and calculate celltype mean for each treatment group
    meta <- as.data.frame(SummarizedExperiment::colData(sce_subset_tmp))
    if (!is.null(splitby_key)) {
        sce_list <- list()
        sce_list_alt <- list()
        for (x in unique(meta[, splitby_key])) {
            sce_list[[x]] <- list()
            sce_list_alt[[x]] <- list()
        }
        for (n in names(sce_list)) {
            for (x in unique(meta[, celltype_key])) {
                sce_list[[n]][[x]] <- sce_subset_tmp[, meta[, celltype_key] == x &
                    meta[, splitby_key] == n]
                sce_list_alt[[n]][[x]] <- sce_subset[, meta[, celltype_key] == x &
                    meta[, splitby_key] == n]
            }
        }
        sce_list2 <- lapply(sce_list, function(y) {
            z <- lapply(y, .cellTypeMeans)
            return(z)
        })
        sce_list3 <- lapply(sce_list, function(y) {
            z <- lapply(y, .cellTypeFraction)
            return(z)
        })
        sce_list2 <- lapply(sce_list2, function(x) do.call(cbind, x))
        sce_list3 <- lapply(sce_list3, function(x) do.call(cbind, x))
        for (n in names(sce_list2)) {
            colnames(sce_list2[[n]]) <- paste0(paste0(n, "_"), colnames(sce_list2[[n]]))
            colnames(sce_list3[[n]]) <- paste0(paste0(n, "_"), colnames(sce_list3[[n]]))
        }
        sce_list2 <- do.call(cbind, sce_list2)
        sce_list3 <- do.call(cbind, sce_list3)
    } else {
        sce_list <- list()
        sce_list_alt <- list()
        for (x in unique(meta[, celltype_key])) {
            sce_list[[x]] <- sce_subset_tmp[, meta[, celltype_key] == x]
            sce_list_alt[[x]] <- sce_subset[, meta[, celltype_key] == x]
        }
        sce_list2 <- lapply(sce_list, .cellTypeMeans)
        sce_list3 <- lapply(sce_list, .cellTypeFraction)
        sce_list2 <- do.call(cbind, sce_list2)
        sce_list3 <- do.call(cbind, sce_list3)
    }
    keep_a <- which(dictionary$gene_a != "")
    keep_b <- which(dictionary$gene_b != "")
    id_a_dict <- dictionary$id_a[keep_a]
    names(id_a_dict) <- dictionary$gene_a[keep_a]
    id_b_dict <- dictionary$id_b[keep_b]
    names(id_b_dict) <- dictionary$gene_b[keep_b]
    id_dict <- c(id_a_dict, id_b_dict)
    id_dict <- id_dict[!duplicated(names(id_dict))]
    humanreadablename <- c()
    for (i in row.names(sce_list2)) {
        humanreadablename <- c(humanreadablename, as.character(unlist(id_dict[i])))
    }
    rownames(sce_list2) <- humanreadablename
    rownames(sce_list3) <- humanreadablename
    decon_subset <- deconvoluted[deconvoluted$complex_name %in% .findComplex(interactions_subset), ]
    if (nrow(decon_subset) > 0) {
        # although multiple rows are returned, really it's the same value for
        # the same complex
        decon_subset <- split(decon_subset, decon_subset$complex_name)
        decon_subset_expr <- lapply(decon_subset, function(x) {
            x <- x[, colnames(x) %in% colnames(sce_list2)]
            x <- colMeans(x)
            return(x)
        })
        if (!is.null(splitby_key)) {
            decon_subset_fraction <- lapply(decon_subset, function(x) {
                z <- unique(x$gene_name)
                test <- lapply(sce_list_alt, function(y) {
                    return(lapply(y, .cellTypeFraction_complex, genes = z, gene_symbol_mapping = gene_symbol_mapping))
                })
                return(test)
            })
            decon_subset_fraction <- lapply(decon_subset_fraction, function(x) {
                y <- lapply(x, function(z) do.call(cbind, z))
                for (i in 1:length(y)) {
                    colnames(y[[i]]) <- paste0(names(y[i]), "_", colnames(y[[i]]))
                }
                y <- do.call(cbind, y)
                return(y)
            })
        } else {
            decon_subset_fraction <- lapply(decon_subset, function(x) {
                z <- unique(x$gene_name)
                test <- lapply(sce_list_alt, function(y) {
                    return(.cellTypeFraction_complex(y, genes = z, gene_symbol_mapping = gene_symbol_mapping))
                })
                return(do.call(cbind, test))
            })
        }
        decon_subset_expr <- do.call(rbind, decon_subset_expr)
        decon_subset_fraction <- do.call(rbind, decon_subset_fraction)
        row.names(decon_subset_fraction) <- row.names(decon_subset_expr)
        requireNamespace("plyr")
        expr_df <- plyr::rbind.fill.matrix(sce_list2, decon_subset_expr)
        row.names(expr_df) <- c(row.names(sce_list2), row.names(decon_subset_expr))
        fraction_df <- plyr::rbind.fill.matrix(sce_list3, decon_subset_fraction)
        row.names(fraction_df) <- c(row.names(sce_list3), row.names(decon_subset_fraction))
    } else {
        expr_df <- sce_list2
        fraction_df <- sce_list3
    }
    # make a big fat edgelist
    if (!is.null(desiredInteractions)) {
        if (class(desiredInteractions) == "list") {
            desiredInteractions_ <- c(desiredInteractions, lapply(
                desiredInteractions,
                rev
            ))
            cell_type_grid <- as.data.frame(do.call(rbind, desiredInteractions_))
        } else if ((class(desiredInteractions) == "data.frame")) {
            cell_type_grid <- desiredInteractions
        }
        cells_test <- unique(unlist(desiredInteractions))
    } else {
        cells_test <- tryCatch(unique(droplevels(meta[, celltype_key])), error = function(e) {
            unique(meta[, celltype_key])
        })
        cell_type_grid <- expand.grid(cells_test, cells_test)
    }
    if (remove_self) {
        rm_idx <- which(cell_type_grid[, 1] == cell_type_grid[, 2])
        if (length(rm_idx) > 0) {
            cell_type_grid <- cell_type_grid[-rm_idx, ]
        }
    }
    ligand <- interactions_subset$id_a
    receptor <- interactions_subset$id_b
    pair <- interactions_subset$interacting_pair
    converted_pair <- interactions_subset$converted
    receptor_a <- interactions_subset$receptor_a
    receptor_b <- interactions_subset$receptor_b
    producers <- as.character(cell_type_grid[, 1])
    receivers <- as.character(cell_type_grid[, 2])
    barcodes <- paste0(lr_interactions$Var2, DEFAULT_SEP, lr_interactions$Var1)
    dfx <- list()
    if (!is.null(splitby_key)) {
        for (i in unique(meta[, splitby_key])) {
            dfx[[i]] <- .generateDf(
                ligand = ligand, sep = DEFAULT_SEP, receptor = receptor,
                receptor_a = receptor_a, receptor_b = receptor_b, pair = pair, converted_pair = converted_pair,
                producers = producers, receivers = receivers, cell_type_means = expr_df,
                cell_type_fractions = fraction_df, sce = sce_subset, sce_alt = sce_list_alt,
                gsm = gene_symbol_mapping, splitted = i
            )
            dfx[[i]] <- dfx[[i]][dfx[[i]]$barcode %in% barcodes, ]
        }
    } else {
        dfx[[1]] <- .generateDf(
            ligand = ligand, sep = DEFAULT_SEP, receptor = receptor,
            receptor_a = receptor_a, receptor_b = receptor_b, pair = pair, converted_pair = converted_pair,
            producers = producers, receivers = receivers, cell_type_means = expr_df,
            cell_type_fractions = fraction_df, sce = sce_subset, sce_alt = sce_list_alt,
            gsm = gene_symbol_mapping
        )
        dfx[[1]] <- dfx[[1]][dfx[[1]]$barcode %in% barcodes, ]
    }

    gl <- list()
    if (length(show_legend) > 1) {
        for (i in 1:length(dfx)) {
            gl[[i]] <- tryCatch(
                .chord_diagram3(
                    tmp_df = dfx[[i]], lr_interaction = lr_interactions,
                    scaled = standard_scale, sep = DEFAULT_SEP, alpha = alpha, directional = directional,
                    show_legend = show_legend[i], edge_cols = edge_colors, grid_cols = grid_colors,
                    legend.pos.x = legend.pos.x, legend.pos.y = legend.pos.y, title = names(dfx)[i]
                ),
                error = function(e) {
                    return(NA)
                }
            )
        }
    } else {
        for (i in 1:length(dfx)) {
            gl[[i]] <- tryCatch(
                .chord_diagram3(
                    tmp_dfx = dfx[[i]], lr_interaction = lr_interactions,
                    scaled = standard_scale, sep = DEFAULT_SEP, alpha = alpha, directional = directional,
                    show_legend = show_legend, edge_cols = edge_colors, grid_cols = grid_colors,
                    legend.pos.x = legend.pos.x, legend.pos.y = legend.pos.y, title = names(dfx)[i]
                ),
                error = function(e) {
                    return(NA)
                }
            )
        }
    }
    if (length(gl) > 1) {
        return(gl)
    } else {
        return(gl[[1]])
    }
}
