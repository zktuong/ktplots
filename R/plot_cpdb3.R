#' Plotting cellphonedb results as a chord diagram
#'
#' @param cell_type1 cell type 1
#' @param cell_type2 cell type 2
#' @param scdata single-cell data. can be seurat/summarizedexperiment object
#' @param idents vector holding the idents for each cell or column name of scdata's metadata. MUST match cpdb's columns
#' @param means object holding means.txt from cpdb output
#' @param pvals object holding pvals.txt from cpdb output
#' @param deconvoluted object holding deconvoluted.txt from cpdb output
#' @param p.adjust.method correction method. p.adjust.methods of one of these options: c('holm', 'hochberg', 'hommel', 'bonferroni', 'BH', 'BY', 'fdr', 'none')
#' @param keep_significant_only logical. Default is FALSE. Switch to TRUE if you only want to plot the significant hits from cpdb.
#' @param split.by column name in the metadata/coldata table to split the spots by. Can only take columns with binary options. If specified, name to split by MUST be specified in the meta file provided to cpdb prior to analysis.
#' @param scale logical. scale the expression to mean +/- SD. NULL defaults to TRUE.
#' @param standard_scale logical. scale the expression to range from 0 to 1. Default is TRUE
#' @param separator default = NULL. separator to use to split between celltypes. Unless otherwise specified, the separator will be `>@<`. Make sure the idents and split.by doesn't overlap with this.
#' @param gene_symbol_mapping default = NULL.column name for rowData in sce holding the actual gene symbols if row names aren't gene symbols
#' @param frac default = 0.1. Minimum fraction of celtypes expressing a gene in order to keep the interaction. Gene must be expressesd >= `frac` in either of the pair of celltypes in order to keep.
#' @param remove_self default = TRUE. Remove self-self arcs.
#' @param desiredInteractions default = NULL. Specific list of celltype comparisons e.g. list(c('CD4_Tcm', 'cDC1'), c('CD4_Tcm', 'cDC2')). Also accepts a dataframe where first column is celltype 1 and 2nd column is celltype 2.
#' @param version3 if is cellphonedb version3. 
#' @param directional Whether links have directions. 1 means the direction is from the first column in df to the second column, -1 is the reverse, 0 is no direction, and 2 for two directional.
#' @param alpha transparency for links
#' @param edge_colors vector of colors for links
#' @param grid_colors vector of colrs for grids
#' @param show_legend whether or not to show the legend
#' @param legend.pos.x x position of legend
#' @param legend.pos.y y position of legend
#' @param ... passes arguments plot_cpdb
#' @return Plotting cellphonedb results as a CellChat inspired chord diagram
#' @examples
#' \donttest{
#' }
#' @importFrom circlize circos.clear chordDiagram
#' @importFrom grDevices recordPlot
#' @export

plot_cpdb3 <- function(cell_type1, cell_type2, scdata, idents, means, pvals, deconvoluted,
    p.adjust.method = NULL, keep_significant_only = TRUE, split.by = NULL, standard_scale = TRUE,
    separator = NULL, gene_symbol_mapping = NULL, frac = 0.1, remove_self = TRUE,
    desiredInteractions = NULL, version3 = FALSE, directional = 1,
    alpha = 0.5, edge_colors = NULL, grid_colors = NULL, show_legend = TRUE, legend.pos.x = 20,
    legend.pos.y = 20, ...) {
    if (class(scdata) == "Seurat") {
        stop("Sorry not supported. Please use a SingleCellExperiment object.")
    }
    if (length(separator) > 0) {
        sep = separator
    } else {
        sep = ">@<"
    }

    lr_interactions = plot_cpdb(cell_type1 = cell_type1, cell_type2 = cell_type2,
        scdata = scdata, idents = idents, split.by = split.by, means = means, pvals = pvals,
        keep_significant_only = keep_significant_only, standard_scale = standard_scale,
        return_table = TRUE, version3 = version3, ...)
    requireNamespace("SummarizedExperiment")
    requireNamespace("SingleCellExperiment")
    if (is.null(split.by)) {
        if (any(lr_interactions[, 3] > 0)) {
            if (any(is.na(lr_interactions[, 3]))) {
                lr_interactions <- lr_interactions[lr_interactions[, 3] > 0 & !is.na(lr_interactions[,
                  3]), ]
            } else {
                lr_interactions <- lr_interactions[lr_interactions[, 3] > 0, ]
            }
        }
    }
    subset_clusters <- unique(unlist(lapply(as.character(lr_interactions$group),
        strsplit, sep)))
    sce_subset <- scdata[, SummarizedExperiment::colData(scdata)[, idents] %in% subset_clusters]
    interactions <- means[, c("interacting_pair", "gene_a", "gene_b", "partner_a",
        "partner_b")]
    interactions$converted <- gsub("-", " ", interactions$interacting_pair)
    interactions$converted <- gsub("_", "-", interactions$converted)
    interactions_subset <- interactions[interactions$converted %in% lr_interactions$Var1,
        ]
    tm0 <- do.call(c, lapply(as.list(interactions_subset$interacting_pair), strsplit,
        "_"))
    if (any(lapply(tm0, length) > 2)) {
        complex_id <- which(lapply(tm0, length) > 2)
        interactions_subset_ = interactions_subset[complex_id, ]
        simple_1 = interactions_subset_$interacting_pair[grep("complex:", interactions_subset_$partner_b)]
        partner_1 = gsub("complex:", "", interactions_subset_$partner_b[grep("complex:",
            interactions_subset_$partner_b)])
        partner_2 = gsub("complex:", "", interactions_subset_$partner_a[grep("complex:",
            interactions_subset_$partner_a)])
        simple_2 = interactions_subset_$interacting_pair[grep("complex:", interactions_subset_$partner_a)]
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
        tm0 = data.frame(t(matrix(unlist(tm0), 2, length(unlist(tm0))/2)))
        colnames(tm0) <- c("id_a", "id_b")
        interactions_subset <- cbind(interactions_subset, tm0)
        dictionary <- interactions_subset[, c("gene_a", "gene_b", "partner_a", "partner_b",
            "id_a", "id_b")]
    } else {
        tm0 = data.frame(t(matrix(unlist(tm0), 2, length(unlist(tm0))/2)))
        colnames(tm0) <- c("id_a", "id_b")
        interactions_subset <- cbind(interactions_subset, tm0)
        dictionary <- interactions_subset[, c("gene_a", "gene_b", "partner_a", "partner_b",
            "id_a", "id_b")]
    }

    # extract all the possible genes.
    geneid = unique(c(interactions_subset$id_a, interactions_subset$id_b))
    # rmg = which(geneid == '') if (length(rmg) > 0){ geneid =
    # geneid[-which(geneid == '')] }
    if (all(!geneid %in% row.names(sce_subset))) {
        geneid = unique(c(interactions_subset$gene_a, interactions_subset$gene_b))
    }
    sce_subset_tmp <- sce_subset[row.names(sce_subset) %in% geneid, ]

    if (dim(sce_subset_tmp)[1] == 0) {
        stop("Gene ids not found. Are you sure your single-cell object contains human gene symbols?")
    }
    # split to list and calculate celltype mean for each treatment group
    meta <- as.data.frame(SummarizedExperiment::colData(sce_subset_tmp))

    cellTypeMeans <- function(x) {
        cm <- Matrix::rowMeans(SingleCellExperiment::counts(x))
        return(cm)
    }

    cellTypeFraction <- function(x) {
        cm <- Matrix::rowMeans(SingleCellExperiment::counts(x) > 0)
        return(cm)
    }

    if (!is.null(split.by)) {
        sce_list <- list()
        sce_list_alt <- list()
        for (x in unique(meta[, split.by])) {
            sce_list[[x]] <- list()
            sce_list_alt[[x]] <- list()
        }

        for (n in names(sce_list)) {
            for (x in unique(meta[, idents])) {
                sce_list[[n]][[x]] <- sce_subset_tmp[, meta[, idents] == x & meta[,
                  split.by] == n]
                sce_list_alt[[n]][[x]] <- sce_subset[, meta[, idents] == x & meta[,
                  split.by] == n]
            }
        }

        sce_list2 <- lapply(sce_list, function(y) {
            z <- lapply(y, cellTypeMeans)
            return(z)
        })
        sce_list3 <- lapply(sce_list, function(y) {
            z <- lapply(y, cellTypeFraction)
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

        for (x in unique(meta[, idents])) {
            sce_list[[x]] <- sce_subset_tmp[, meta[, idents] == x]
            sce_list_alt[[x]] <- sce_subset[, meta[, idents] == x]
        }

        sce_list2 <- lapply(sce_list, cellTypeMeans)
        sce_list3 <- lapply(sce_list, cellTypeFraction)

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
    humanreadablename = c()
    for (i in row.names(sce_list2)) {
        humanreadablename = c(humanreadablename, as.character(unlist(id_dict[i])))
    }
    rownames(sce_list2) <- humanreadablename
    rownames(sce_list3) <- humanreadablename

    findComplex <- function(interaction) {
        idxa <- which(interaction$gene_a == "")
        idxb <- which(interaction$gene_b == "")
        complexa <- gsub("complex:", "", interaction$partner_a[idxa])
        complexb <- gsub("complex:", "", interaction$partner_b[idxb])

        if (length(complexa) > 0) {
            if (length(complexb) > 0) {
                res <- c(complexa, complexb)
            } else {
                res <- complexa
            }
        } else if (length(complexb) > 0) {
            res <- complexb
        } else {
            res <- NULL
        }
        return(res)
    }

    # Utility function to retrieve the mean
    cellTypeExpr_complex <- function(sce_, genes, gene_symbol_mapping = NULL) {
        scex <- tryCatch(sce_[genes, ], error = function(e) {
            if (!is.null(gene_symbol_mapping)) {
                sce_[which(rowData(sce_)[, gene_symbol_mapping] %in% genes), ]
            } else {
                sce_[which(rowData(sce_)[, "index"] %in% genes), ]
            }
        })

        cm <- mean(Matrix::rowMeans(SingleCellExperiment::counts(scex)))
        return(cm)
    }
    # to retrieve the fraction, we use the average fraction of all genes
    # mapping to the complex
    cellTypeFraction_complex <- function(sce_, genes, gene_symbol_mapping = NULL) {
        scex <- tryCatch(sce_[genes, ], error = function(e) {
            if (!is.null(gene_symbol_mapping)) {
                sce_[which(rowData(sce_)[, gene_symbol_mapping] %in% genes), ]
            } else {
                sce_[which(rowData(sce_)[, "index"] %in% genes), ]
            }
        })

        cm <- mean(Matrix::rowMeans(SingleCellExperiment::counts(scex) > 0))
        return(cm)
    }

    decon_subset <- deconvoluted[deconvoluted$complex_name %in% findComplex(interactions_subset),
        ]
    if (nrow(decon_subset) > 0) {
        # although multiple rows are returned, really it's the same value for
        # the same complex
        decon_subset <- split(decon_subset, decon_subset$complex_name)
        decon_subset_expr <- lapply(decon_subset, function(x) {
            x <- x[, colnames(x) %in% colnames(sce_list2)]
            x <- colMeans(x)
            return(x)
        })

        if (!is.null(split.by)) {
            decon_subset_fraction <- lapply(decon_subset, function(x) {
                x <- unique(x$gene_name)
                test <- lapply(sce_list_alt, function(y) {
                  return(lapply(y, cellTypeFraction_complex, x, gene_symbol_mapping))
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
                x <- unique(x$gene_name)
                test <- lapply(sce_list_alt, function(y) {
                  return(cellTypeFraction_complex(y, x, gene_symbol_mapping))
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
            desiredInteractions_ <- c(desiredInteractions, lapply(desiredInteractions,
                rev))
            cell_type_grid <- as.data.frame(do.call(rbind, desiredInteractions_))
        } else if ((class(desiredInteractions) == "data.frame")) {
            cell_type_grid <- desiredInteractions
        }
        cells_test = unique(unlist(desiredInteractions))
    } else {
        cells_test <- tryCatch(unique(droplevels(meta[, idents])), error = function(e) unique(meta[,
            idents]))
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
    producers <- as.character(cell_type_grid[, 1])
    receivers <- as.character(cell_type_grid[, 2])

    generateDf <- function(ligand, sep, receptor, pair, converted_pair, producers,
        receivers, cell_type_means, cell_type_fractions, splitted = NULL) {
        if (!is.null(splitted)) {
            pp <- paste0(splitted, "_", producers)
            rc <- paste0(splitted, "_", receivers)
        } else {
            pp <- producers
            rc <- receivers
        }
        producer_expression <- data.frame()
        producer_fraction <- data.frame()
        for (i in seq_along(pp)) {
            for (j in seq_along(ligand)) {
                if (any(grepl(paste0("^", ligand[j], "$"), row.names(cell_type_means)))) {
                  x <- cell_type_means[ligand[j], pp[i]]
                  y <- cell_type_fractions[ligand[j], pp[i]]
                } else {
                  if (any(grepl(paste0("^", ligand[j], "$"), row.names(sce_subset)))) {
                    x <- cellTypeExpr_complex(sce_list_alt[[pp[i]]], ligand[j], gene_symbol_mapping)
                    y <- cellTypeFraction_complex(sce_list_alt[[pp[i]]], ligand[j],
                      gene_symbol_mapping)
                  } else {
                    x <- 0
                    y <- 0
                  }
                }
                producer_expression[ligand[j], pp[i]] <- x
                producer_fraction[ligand[j], pp[i]] <- y
            }
        }

        receiver_expression <- data.frame()
        receiver_fraction <- data.frame()
        for (i in seq_along(rc)) {
            for (j in seq_along(receptor)) {
                if (any(grepl(paste0("^", receptor[j], "$"), row.names(cell_type_means)))) {
                  x <- cell_type_means[receptor[j], rc[i]]
                  y <- cell_type_fractions[receptor[j], rc[i]]
                } else {
                  if (any(grepl(paste0("^", receptor[j], "$"), row.names(sce_subset)))) {
                    x <- cellTypeExpr_complex(sce_list_alt[[rc[i]]], receptor[j],
                      gene_symbol_mapping)
                    y <- cellTypeFraction_complex(sce_list_alt[[rc[i]]], receptor[j],
                      gene_symbol_mapping)
                  } else {
                    x <- 0
                    y <- 0
                  }
                }
                receiver_expression[receptor[j], rc[i]] <- x
                receiver_fraction[receptor[j], rc[i]] <- y
            }
        }
        test_df = list()
        for (i in seq_along(pp)) {
            px <- pp[i]
            rx <- rc[i]
            for (j in seq_along(pair)) {
                lg <- ligand[j]
                rcp <- receptor[j]
                pr <- pair[j]
                out = data.frame(c(lg, rcp, pr, px, rx, producer_expression[lg, px],
                  producer_fraction[lg, px], receiver_expression[rcp, rx], receiver_fraction[rcp,
                    rx]))
                test_df = c(test_df, out)
            }
        }
        df_ <- do.call(rbind, test_df)
        row.names(df_) <- 1:nrow(df_)
        colnames(df_) <- c("ligand", "receptor", "pair", "producer", "receiver",
            "producer_expression", "producer_fraction", "receiver_expression", "receiver_fraction")
        df_ <- as.data.frame(df_)
        df_$from = paste0(df_$producer, sep, df_$ligand)
        df_$to = paste0(df_$receiver, sep, df_$receptor)
        if (!is.null(splitted)) {
            df_$producer_ = df_$producer
            df_$receiver_ = df_$receiver
            df_$from = gsub(paste0(splitted, "_"), "", df_$from)
            df_$to = gsub(paste0(splitted, "_"), "", df_$to)
            df_$producer = gsub(paste0(splitted, "_"), "", df_$producer)
            df_$receiver = gsub(paste0(splitted, "_"), "", df_$receiver)
            df_$barcode = paste0(df_$producer_, "-", df_$receiver_, sep, converted_pair)
        } else {
            df_$barcode = paste0(df_$producer, "-", df_$receiver, sep, converted_pair)
        }

        return(df_)
    }
    barcodes = paste0(lr_interactions$Var2, sep, lr_interactions$Var1)
    dfx <- list()
    if (!is.null(split.by)) {
        for (i in unique(meta[, split.by])) {
            dfx[[i]] <- generateDf(ligand, sep, receptor, pair, converted_pair, producers,
                receivers, expr_df, fraction_df, i)
            dfx[[i]] <- dfx[[i]][dfx[[i]]$barcode %in% barcodes, ]
        }
    } else {
        dfx[[1]] = generateDf(ligand, sep, receptor, pair, converted_pair, producers,
            receivers, expr_df, fraction_df)
        dfx[[1]] <- dfx[[1]][dfx[[1]]$barcode %in% barcodes, ]
    }

    scPalette <- function(n) {
        colorSpace <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#F29403", "#F781BF",
            "#BC9DCC", "#A65628", "#54B0E4", "#222F75", "#1B9E77", "#B2DF8A", "#E3BE00",
            "#FB9A99", "#E7298A", "#910241", "#00CDD1", "#A6CEE3", "#CE1261", "#5E4FA2",
            "#8CA77B", "#00441B", "#DEDC00", "#B3DE69", "#8DD3C7", "#999999")
        if (n <= length(colorSpace)) {
            colors <- colorSpace[1:n]
        } else {
            colors <- (grDevices::colorRampPalette(colorSpace))(n)
        }
        return(colors)
    }

    chord_diagram <- function(tmp_dfx, lr_interactions, p.adjust_method, scaled,
        alpha, directional, show_legend, edge_cols, grid_cols, legend.pos.x, legend.pos.y) {

        if (scaled) {
            interactions_items <- lr_interactions$scaled_means
        } else {
            interactions_items <- lr_interactions$means
        }
        names(interactions_items) <- paste0(lr_interactions$Var2, sep, lr_interactions$Var1)
        if (!is.null(p.adjust_method)) {
            pvals_items <- lr_interactions$pvals_adj
        } else {
            pvals_items <- lr_interactions$pvals
        }
        names(pvals_items) <- paste0(lr_interactions$Var2, sep, lr_interactions$Var1)
        interactions_items[is.na(pvals_items)] <- 1
        tmp_dfx$pair <- gsub("_", " - ", tmp_dfx$pair)
        tmp_dfx$value <- interactions_items[tmp_dfx$barcode]
        tmp_dfx$pval <- pvals_items[tmp_dfx$barcode]
        if (!is.null(edge_cols)) {
            if (length(edge_cols) != length(unique(tmp_dfx$pair))) {
                stop(paste0("Please provide ", length(unique(tmp_dfx$pair)), " to edge_colors."))
            } else {
                edge_color <- edge_cols
            }
        } else {
            edge_color <- scPalette(length(unique(tmp_dfx$pair)))
        }

        if (!is.null(grid_cols)) {
            if (length(grid_cols) != length(unique(tmp_dfx$receiver))) {
                stop(paste0("Please provide ", length(unique(tmp_dfx$receiver)),
                  " to grid_colors."))
            } else {
                grid_color <- grid_cols
            }
        } else {
            grid_color <- scPalette(length(unique(tmp_dfx$receiver)))
        }
        names(edge_color) <- unique(tmp_dfx$pair)
        names(grid_color) <- unique(tmp_dfx$receiver)
        tmp_dfx$edge_color = edge_color[tmp_dfx$pair]
        tmp_dfx$edge_color <- colorspace::adjust_transparency(tmp_dfx$edge_color,
            alpha = alpha)
        tmp_dfx$edge_color[is.na(tmp_dfx$pval)] <- "#00000000"
        tmp_dfx$grid_color = grid_color[tmp_dfx$receiver]

        tmp_dfx <- tmp_dfx[!duplicated(tmp_dfx$barcode), ]

        if (directional == 2) {
            link.arr.type = "triangle"
        } else {
            link.arr.type = "big.arrow"
        }

        cells <- unique(c(tmp_dfx$producer, tmp_dfx$receiver))
        names(cells) <- cells

        circos.clear()
        chordDiagram(tmp_dfx[c("producer", "receiver", "value")], directional = directional,
            direction.type = c("diffHeight", "arrows"), link.arr.type = link.arr.type,
            annotationTrack = c("name", "grid"), col = tmp_dfx$edge_color, grid.col = grid_color,
            group = cells)
        requireNamespace("ComplexHeatmap")
        if (show_legend) {
            lgd <- ComplexHeatmap::Legend(at = names(edge_color), type = "grid",
                legend_gp = grid::gpar(fill = edge_color), title = "interactions")
            ComplexHeatmap::draw(lgd, x = unit(1, "npc") - unit(legend.pos.x, "mm"),
                y = unit(legend.pos.y, "mm"), just = c("right", "bottom"))
        }

        circos.clear()

        gg <- recordPlot()
        return(gg)
    }
    gl <- list()
    if (length(show_legend) > 1) {
        for (i in 1:length(dfx)) {
            gl[[i]] <- tryCatch(chord_diagram(dfx[[i]], lr_interactions, p.adjust.method,
                standard_scale, alpha, directional, show_legend[i], edge_colors,
                grid_colors, legend.pos.x, legend.pos.y), error = function(e) return(NA))
        }
    } else {
        for (i in 1:length(dfx)) {
            gl[[i]] <- tryCatch(chord_diagram(dfx[[i]], lr_interactions, p.adjust.method,
                standard_scale, alpha, directional, show_legend, edge_colors, grid_colors,
                legend.pos.x, legend.pos.y), error = function(e) return(NA))
        }
    }
    
    if (length(gl) > 1) {
        return(gl)
    } else {
        return(gl[[1]])
    }

}