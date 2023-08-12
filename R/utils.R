#' @import ggplot2
#' @import ggraph
#' @importFrom circlize circos.clear chordDiagram
#' @importFrom grDevices recordPlot

DEFAULT_SEP <- ">@<"
DEFAULT_SPEC_PAT <- "/|:|\\?|\\*|\\+|[\\]|\\(|\\)|\\/"


.prep_table <- function(data) {
    dat <- data
    rownames(dat) <- make.names(dat$interacting_pair, unique = TRUE)
    colnames(dat) <- gsub("\\|", DEFAULT_SEP, colnames(dat))
    rownames(dat) <- gsub("_", "-", rownames(dat))
    rownames(dat) <- gsub("[.]", " ", rownames(dat))
    return(dat)
}

.set_x_stroke <- function(df, isnull, stroke) {
    for (i in seq_len(nrow(df))) {
        if (isnull) {
            nullstatus <- is.na(df[i, "x_stroke"])
        } else {
            nullstatus <- !is.na(df[i, "x_stroke"])
        }
        if (nullstatus) {
            df[i, "x_stroke"] <- stroke
        }
    }
    return(df)
}

.prep_query_group <- function(data, genes = NULL, gene_family = NULL, custom_gene_family = NULL) {
    if (is.null(gene_family) & is.null(genes)) {
        query_group <- NULL
        query <- grep("", data$interacting_pair)
    } else if (!is.null(gene_family) & !is.null(genes)) {
        stop("Please specify either genes or gene_family, not both")
    } else if (!is.null(gene_family) & is.null(genes)) {
        chemokines <- grep("^CXC|CCL|CCR|CX3|XCL|XCR", data$interacting_pair)
        th1 <- grep("IL2|IL12|IL18|IL27|IFNG|IL10|TNF$|TNF |LTA|LTB|STAT1|CCR5|CXCR3|IL12RB1|IFNGR1|TBX21|STAT4", data$interacting_pair)
        th2 <- grep("IL4|IL5|IL25|IL10|IL13|AREG|STAT6|GATA3|IL4R", data$interacting_pair)
        th17 <- grep("IL21|IL22|IL24|IL26|IL17A|IL17A|IL17F|IL17RA|IL10|RORC|RORA|STAT3|CCR4|CCR6|IL23RA|TGFB", data$interacting_pair)
        treg <- grep("IL35|IL10|FOXP3|IL2RA|TGFB", data$interacting_pair)
        costimulatory <- grep("CD86|CD80|CD48|LILRB2|LILRB4|TNF|CD2|ICAM|SLAM|LT[AB]|NECTIN2|CD40|CD70|CD27|CD28|CD58|TSLP|PVR|CD44|CD55|CD[1-9]", data$interacting_pair)
        coinhibitory <- grep("SIRP|CD47|ICOS|TIGIT|CTLA4|PDCD1|CD274|LAG3|HAVCR|VSIR", data$interacting_pair)
        query_group <- list(
            chemokines = chemokines,
            chemokine = chemokines,
            th1 = th1,
            th2 = th2,
            th17 = th17,
            treg = treg,
            costimulatory = costimulatory,
            coinhibitory = coinhibitory,
            costimulation = costimulatory,
            coinhibition = coinhibitory
        )
        if (!is.null(custom_gene_family)) {
            cgf <- as.list(custom_gene_family)
            cgf <- lapply(cgf, function(x) grep(paste(x, collapse = "|"), data$interacting_pair))
            query_group <- c(query_group, cgf)
        }
        query <- NULL
    } else if (is.null(gene_family) & !is.null(genes)) {
        query_group <- NULL
        query <- grep(paste(genes, collapse = "|"), data$interacting_pair)
    }
    out <- list("query_group" = query_group, "query" = query)
    return(out)
}

.sub_pattern <- function(cell_type, pattern) {
    cell_type_tmp <- unlist(strsplit(cell_type, "*"))
    if (any(grepl(pattern, cell_type_tmp))) {
        idxz <- grep(pattern, cell_type_tmp)
        cell_type_tmp[idxz] <- paste0("\\", cell_type_tmp[idxz])
        cell_typex <- do.call(paste, c(as.list(cell_type_tmp), sep = ""))
    } else {
        cell_typex <- cell_type
    }
    return(cell_typex)
}

.gg_color_hue <- function(n) {
    requireNamespace("grDevices")
    hues <- seq(15, 375, length = n + 1)
    grDevices::hcl(h = hues, l = 65, c = 100)[1:n]
}

.prep_data_querygroup_celltype1 <- function(.data, .query_group, .gene_family, .cell_type, .celltype, ...) {
    dat <- suppressWarnings(tryCatch(
        .data[.query_group[[tolower(.gene_family)]],
            grep(.cell_type, colnames(.data), useBytes = TRUE, ...),
            drop = FALSE
        ],
        error = function(e) {
            colidx <- lapply(.celltype, function(z) {
                grep(z, colnames(.data),
                    useBytes = TRUE, ...
                )
            })
            colidx <- unique(do.call(c, colidx))
            tmpm <- .data[.query_group[[tolower(.gene_family)]], colidx, drop = FALSE]
            return(tmpm)
        }
    ))
    return(dat)
}

.prep_data_querygroup_celltype2 <- function(.data, .query_group, .gene_family, .cell_type, .celltype, ...) {
    dat <- suppressWarnings(tryCatch(
        .data[unlist(.query_group[c(tolower(.gene_family))], use.names = FALSE),
            grep(.cell_type, colnames(.data), useBytes = TRUE, ...),
            drop = FALSE
        ],
        error = function(e) {
            colidx <- lapply(.celltype, function(z) {
                grep(z, colnames(.data),
                    useBytes = TRUE, ...
                )
            })
            colidx <- unique(do.call(c, colidx))
            tmpm <- .data[unlist(.query_group[c(tolower(.gene_family))], use.names = FALSE), colidx, drop = FALSE]
            return(tmpm)
        }
    ))
    return(dat)
}

.prep_data_query_celltype <- function(.data, .query, .cell_type, .celltype, ...) {
    dat <- suppressWarnings(tryCatch(.data[.query, grep(.cell_type, colnames(.data),
        useBytes = TRUE, ...
    ), drop = FALSE], error = function(e) {
        colidx <- lapply(.celltype, function(z) {
            grep(z, colnames(.data),
                useBytes = TRUE,
                ...
            )
        })
        colidx <- unique(do.call(c, colidx))
        tmpm <- .data[.query, colidx, drop = FALSE]
        return(tmpm)
    }))
    return(dat)
}


.create_celltype_query <- function(ctype1, ctype2, sep) {
    ct1 <- list()
    ct2 <- list()
    for (i in 1:length(ctype2)) {
        ct1[i] <- paste0("^", ctype1, sep, ctype2[i], "$")
        ct2[i] <- paste0("^", ctype2[i], sep, ctype1, "$")
    }
    ct_1 <- do.call(paste0, list(ct1, collapse = "|"))
    ct_2 <- do.call(paste0, list(ct2, collapse = "|"))
    ct <- list(ct_1, ct_2)
    ct <- do.call(paste0, list(ct, collapse = "|"))
    return(ct)
}

.keep_interested_groups <- function(g, ct, sep) {
    ctx <- strsplit(ct, "\\|")[[1]]
    idx <- grep(paste0(g, paste0(".*", sep), g), ctx)
    ctx <- ctx[idx]
    ctx <- paste0(ctx, collapse = "|")
    return(ctx)
}

.scPalette <- function(n) {
    requireNamespace("grDevices")
    colorSpace <- c(
        "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#F29403", "#F781BF",
        "#BC9DCC", "#A65628", "#54B0E4", "#222F75", "#1B9E77", "#B2DF8A", "#E3BE00",
        "#FB9A99", "#E7298A", "#910241", "#00CDD1", "#A6CEE3", "#CE1261", "#5E4FA2",
        "#8CA77B", "#00441B", "#DEDC00", "#B3DE69", "#8DD3C7", "#999999"
    )
    if (n <= length(colorSpace)) {
        colors <- colorSpace[1:n]
    } else {
        colors <- (grDevices::colorRampPalette(colorSpace))(n)
    }
    return(colors)
}

.cellTypeMeans <- function(x) {
    requireNamespace("Matrix")
    requireNamespace("SingleCellExperiment")
    cm <- Matrix::rowMeans(SingleCellExperiment::counts(x))
    return(cm)
}

.cellTypeFraction <- function(x) {
    requireNamespace("Matrix")
    requireNamespace("SingleCellExperiment")
    cm <- Matrix::rowMeans(SingleCellExperiment::counts(x) > 0)
    return(cm)
}

.findComplex <- function(interaction) {
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

# Utility function to retrieve the mean for complex
.cellTypeExpr_complex <- function(sce_, genes, gene_symbol_mapping = NULL) {
    requireNamespace("Matrix")
    requireNamespace("SingleCellExperiment")
    scex <- tryCatch(sce_[genes, ], error = function(e) {
        if (!is.null(gene_symbol_mapping)) {
            sce_[which(SingleCellExperiment::rowData(sce_)[, gene_symbol_mapping] %in%
                genes), ]
        } else {
            sce_[which(SingleCellExperiment::rowData(sce_)[, "index"] %in% genes), ]
        }
    })
    cm <- mean(Matrix::rowMeans(SingleCellExperiment::counts(scex)))
    return(cm)
}

# Utility function to retrieve the fraction for complex
.cellTypeFraction_complex <- function(sce_, genes, gene_symbol_mapping = NULL) {
    requireNamespace("Matrix")
    requireNamespace("SingleCellExperiment")
    scex <- tryCatch(sce_[genes, ], error = function(e) {
        if (!is.null(gene_symbol_mapping)) {
            sce_[which(SingleCellExperiment::rowData(sce_)[, gene_symbol_mapping] %in%
                genes), ]
        } else {
            sce_[which(SingleCellExperiment::rowData(sce_)[, "index"] %in% genes), ]
        }
    })
    cm <- mean(Matrix::rowMeans(SingleCellExperiment::counts(scex) > 0))
    return(cm)
}

.generateDf <- function(
    ligand, sep, receptor, receptor_a, receptor_b, pair, converted_pair,
    producers, receivers, cell_type_means, cell_type_fractions, sce, sce_alt, gsm,
    splitted = NULL) {
    if (!is.null(splitted)) {
        pp <- paste0(splitted, "_", producers)
        rc <- paste0(splitted, "_", receivers)
    } else {
        pp <- producers
        rc <- receivers
    }
    producer_expression <- data.frame()
    producer_fraction <- data.frame()
    if (!is.null(splitted)) {
        sce_altx <- sce_alt[[splitted]]
    } else {
        sce_altx <- sce_alt
    }
    for (i in seq_along(pp)) {
        for (j in seq_along(ligand)) {
            if (any(grepl(paste0("^", ligand[j], "$"), row.names(cell_type_means)))) {
                x <- cell_type_means[ligand[j], pp[i]]
                y <- cell_type_fractions[ligand[j], pp[i]]
            } else {
                if (any(grepl(paste0("^", ligand[j], "$"), row.names(sce)))) {
                    x <- .cellTypeExpr_complex(
                        sce_altx[[producers[i]]],
                        ligand[j], gsm
                    )
                    y <- .cellTypeFraction_complex(
                        sce_altx[[producers[i]]],
                        ligand[j], gsm
                    )
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
                if (any(grepl(paste0("^", receptor[j], "$"), row.names(sce)))) {
                    x <- .cellTypeExpr_complex(
                        sce_altx[[receivers[i]]],
                        receptor[j], gsm
                    )
                    y <- .cellTypeFraction_complex(
                        sce_altx[[receivers[i]]],
                        receptor[j], gsm
                    )
                } else {
                    x <- 0
                    y <- 0
                }
            }
            receiver_expression[receptor[j], rc[i]] <- x
            receiver_fraction[receptor[j], rc[i]] <- y
        }
    }
    test_df <- list()
    for (i in seq_along(pp)) {
        px <- pp[i]
        rx <- rc[i]
        for (j in seq_along(pair)) {
            lg <- ligand[j]
            rcp <- receptor[j]
            ra <- receptor_a[j]
            rb <- receptor_b[j]
            pr <- pair[j]
            out <- data.frame(c(lg, rcp, ra, rb, pr, px, rx, producer_expression[
                lg,
                px
            ], producer_fraction[lg, px], receiver_expression[rcp, rx], receiver_fraction[
                rcp,
                rx
            ]))
            test_df <- c(test_df, out)
        }
    }
    df_ <- do.call(rbind, test_df)
    row.names(df_) <- 1:nrow(df_)
    colnames(df_) <- c(
        "ligand", "receptor", "receptor_a", "receptor_b", "pair",
        "producer", "receiver", "producer_expression", "producer_fraction", "receiver_expression",
        "receiver_fraction"
    )
    df_ <- as.data.frame(df_)
    df_$from <- paste0(df_$producer, sep, df_$ligand)
    df_$to <- paste0(df_$receiver, sep, df_$receptor)
    if (!is.null(splitted)) {
        df_$producer_ <- df_$producer
        df_$receiver_ <- df_$receiver
        df_$from <- gsub(paste0(splitted, "_"), "", df_$from)
        df_$to <- gsub(paste0(splitted, "_"), "", df_$to)
        df_$producer <- gsub(paste0(splitted, "_"), "", df_$producer)
        df_$receiver <- gsub(paste0(splitted, "_"), "", df_$receiver)
        df_$barcode <- paste0(df_$producer_, "-", df_$receiver_, sep, converted_pair)
    } else {
        df_$barcode <- paste0(df_$producer, "-", df_$receiver, sep, converted_pair)
    }
    return(df_)
}
.swap_ligand_receptor <- function(df) {
    is_r_a <- as.logical(df$receptor_a)
    is_r_b <- as.logical(df$receptor_b)
    lg <- df$ligand
    rp <- df$receptor
    from <- df$from
    to <- df$to
    prd <- df$producer
    rec <- df$receiver
    prd_exp <- df$producer_expression
    prd_fra <- df$producer_fraction
    rec_exp <- df$receiver_expression
    rec_fra <- df$receiver_fraction
    # create swaps
    lg_swap <- c()
    rp_swap <- c()
    from_swap <- c()
    to_swap <- c()
    prd_swap <- c()
    rec_swap <- c()
    prd_exp_swap <- c()
    prd_fra_swap <- c()
    rec_exp_swap <- c()
    rec_fra_swap <- c()
    for (i in seq_along(is_r_a)) {
        if (!is_r_a[i]) {
            if (is_r_b[i]) {
                lg_swap <- c(lg_swap, lg[i])
                rp_swap <- c(rp_swap, rp[i])
                from_swap <- c(from_swap, from[i])
                to_swap <- c(to_swap, to[i])
                prd_swap <- c(prd_swap, prd[i])
                rec_swap <- c(rec_swap, rec[i])
                prd_exp_swap <- c(prd_exp_swap, prd_exp[i])
                prd_fra_swap <- c(prd_fra_swap, prd_fra[i])
                rec_exp_swap <- c(rec_exp_swap, rec_exp[i])
                rec_fra_swap <- c(rec_fra_swap, rec_fra[i])
            } else {
                lg_swap <- c(lg_swap, lg[i])
                rp_swap <- c(rp_swap, rp[i])
                from_swap <- c(from_swap, from[i])
                to_swap <- c(to_swap, to[i])
                prd_swap <- c(prd_swap, prd[i])
                rec_swap <- c(rec_swap, rec[i])
                prd_exp_swap <- c(prd_exp_swap, prd_exp[i])
                prd_fra_swap <- c(prd_fra_swap, prd_fra[i])
                rec_exp_swap <- c(rec_exp_swap, rec_exp[i])
                rec_fra_swap <- c(rec_fra_swap, rec_fra[i])
            }
        } else if (is_r_a[i]) {
            if (is_r_b[i]) {
                lg_swap <- c(lg_swap, lg[i])
                rp_swap <- c(rp_swap, rp[i])
                from_swap <- c(from_swap, from[i])
                to_swap <- c(to_swap, to[i])
                prd_swap <- c(prd_swap, prd[i])
                rec_swap <- c(rec_swap, rec[i])
                prd_exp_swap <- c(prd_exp_swap, prd_exp[i])
                prd_fra_swap <- c(prd_fra_swap, prd_fra[i])
                rec_exp_swap <- c(rec_exp_swap, rec_exp[i])
                rec_fra_swap <- c(rec_fra_swap, rec_fra[i])
            } else {
                lg_swap <- c(lg_swap, rp[i])
                rp_swap <- c(rp_swap, lg[i])
                from_swap <- c(from_swap, to[i])
                to_swap <- c(to_swap, from[i])
                prd_swap <- c(prd_swap, rec[i])
                rec_swap <- c(rec_swap, prd[i])
                prd_exp_swap <- c(prd_exp_swap, rec_exp[i])
                prd_fra_swap <- c(prd_fra_swap, rec_fra[i])
                rec_exp_swap <- c(rec_exp_swap, prd_exp[i])
                rec_fra_swap <- c(rec_fra_swap, prd_fra[i])
            }
        }
    }
    df$ligand_swap <- lg_swap
    df$receptor_swap <- rp_swap
    df$pair_swap <- paste0(lg_swap, " - ", rp_swap)
    df$producer_swap <- prd_swap
    df$receiver_swap <- rec_swap
    df$producer_expression_swap <- prd_exp_swap
    df$producer_fraction_swap <- prd_fra_swap
    df$receiever_expression_swap <- rec_exp_swap
    df$receiever_fraction_swap <- rec_fra_swap
    df$from_swap <- from_swap
    df$to_swap <- to_swap
    return(df)
}

.constructGraph <- function(input_group, sep, el, el0, unique_id, interactions_df,
                            plot_cpdb_out, celltype_key, edge_group = FALSE, edge_group_colors = NULL, node_group_colors = NULL, plot_score_as_thickness = TRUE) {
    requireNamespace("igraph")
    celltypes <- unique(c(as.character(el$producer), as.character(el$receiver)))
    el1 <- data.frame(
        from = "root", to = celltypes, barcode_1 = NA, barcode_2 = NA,
        barcode_3 = NA
    )
    el2 <- data.frame(
        from = celltypes, to = paste0(celltypes, sep, "ligand"),
        barcode_1 = NA, barcode_2 = NA, barcode_3 = NA
    )
    el3 <- data.frame(
        from = celltypes, to = paste0(celltypes, sep, "receptor"),
        barcode_1 = NA, barcode_2 = NA, barcode_3 = NA
    )
    el4 <- do.call(rbind, lapply(celltypes, function(x) {
        cell_ligands <- grep(x, el$from, value = TRUE)
        cell_ligands_idx <- grep(x, el$from)
        if (length(cell_ligands) > 0) {
            df <- data.frame(
                from = paste0(x, sep, "ligand"), to = cell_ligands,
                barcode_1 = el$barcode[cell_ligands_idx], barcode_2 = el$pair[cell_ligands_idx],
                barcode_3 = paste0(el$from[cell_ligands_idx], sep, el$to[cell_ligands_idx])
            )
        } else {
            df <- NULL
        }
    }))
    el5 <- do.call(rbind, lapply(celltypes, function(x) {
        cell_ligands <- grep(x, el$to, value = TRUE)
        cell_ligands_idx <- grep(x, el$to)
        if (length(cell_ligands) > 0) {
            df <- data.frame(
                from = paste0(x, sep, "receptor"), to = cell_ligands,
                barcode_1 = el$barcode[cell_ligands_idx], barcode_2 = el$pair[cell_ligands_idx],
                barcode_3 = paste0(el$from[cell_ligands_idx], sep, el$to[cell_ligands_idx])
            )
        } else {
            df <- NULL
        }
    }))
    gr_el <- do.call(rbind, list(el1, el2, el3, el4, el5))
    plot_cpdb_out$barcode <- paste0(plot_cpdb_out$Var2, sep, plot_cpdb_out$Var1)
    mean_col <- grep("means$", colnames(plot_cpdb_out), value = TRUE)
    means <- plot_cpdb_out[
        match(gr_el$barcode_1, plot_cpdb_out$barcode),
        mean_col
    ]
    pval_col <- grep("pvals", colnames(plot_cpdb_out), value = TRUE)
    pvals <- plot_cpdb_out[
        match(gr_el$barcode_1, plot_cpdb_out$barcode),
        pval_col
    ]
    gr_el <- cbind(gr_el, means, pvals)
    if (edge_group) {
        groups <- interactions_df$group[match(gr_el$barcode_2, interactions_df$interacting_pair)]
    }
    gr <- igraph::graph_from_edgelist(as.matrix(gr_el[, 1:2]))
    igraph::E(gr)$interaction_score <- as.numeric(means)
    igraph::E(gr)$pvals <- as.numeric(pvals)
    if (edge_group) {
        igraph::E(gr)$group <- groups
    }
    igraph::E(gr)$name <- gr_el$barcode_3
    # order the graph vertices
    igraph::V(gr)$type <- NA
    igraph::V(gr)$type[igraph::V(gr)$name %in% el4$to] <- "ligand"
    igraph::V(gr)$type[igraph::V(gr)$name %in% el5$to] <- "receptor"
    from <- match(el0$from, igraph::V(gr)$name)
    to <- match(el0$to, igraph::V(gr)$name)
    dat <- data.frame(from = el0$from, to = el0$to)
    if (nrow(dat) > 0) {
        dat$barcode <- paste0(dat$from, sep, dat$to)
        interaction_score <- igraph::E(gr)$interaction_score[match(dat$barcode, gr_el$barcode_3)]
        pval <- igraph::E(gr)$pvals[match(dat$barcode, gr_el$barcode_3)]
        if (any(is.na(pval))) {
            pval[is.na(pval)] <- 1
        }
        if (!all(is.na(range01(-log10(pval))))) {
            pval <- range01(-log10(pval))
        }
        if (edge_group) {
            group <- igraph::E(gr)$group[match(dat$barcode, gr_el$barcode_3)]
        }
        ligand_expr <- data.frame(
            cell_mol = el$from, expression = el$producer_expression,
            fraction = el$producer_fraction
        )
        recep_expr <- data.frame(
            cell_mol = el$to, expression = el$receiver_expression,
            fraction = el$receiver_fraction
        )
        expression <- rbind(ligand_expr, recep_expr)
        df <- igraph::as_data_frame(gr, "both")
        df$vertices$expression <- 0
        df$vertices$fraction <- 0
        df$vertices$expression <- as.numeric(expression$expression)[match(
            df$vertices$name,
            expression$cell_mol
        )]
        df$vertices$fraction <- as.numeric(expression$fraction)[match(
            df$vertices$name,
            expression$cell_mol
        )]
        df$vertices$celltype <- ""
        for (x in unique_id) {
            idx <- grepl(paste0(x, sep), df$vertices$name)
            df$vertices$celltype[idx] <- x
        }
        df$vertices$label <- df$vertices$name
        df$vertices$label[!df$vertices$name %in% c(el0$from, el0$to)] <- ""
        requireNamespace("igraph")
        gr <- igraph::graph_from_data_frame(df$edges, directed = TRUE, vertices = df$vertices)
        for (x in unique_id) {
            igraph::V(gr)$label <- gsub(paste0(x, sep), "", igraph::V(gr)$label)
        }
        if (!is.null(edge_group_colors)) {
            edge_group_colors <- edge_group_colors
        } else {
            nn <- length(unique(igraph::E(gr)$group))
            edge_group_colors <- .gg_color_hue(nn)
        }
        if (!is.null(node_group_colors)) {
            node_group_colors <- node_group_colors
        } else {
            nn <- length(unique(meta[, celltype_key]))
            node_group_colors <- .gg_color_hue(nn)
        }
        # plot the graph
        if (edge_group) {
            if (plot_score_as_thickness) {
                pl <- ggraph(gr, layout = "dendrogram", circular = TRUE) +
                    geom_conn_bundle(
                        data = get_con(
                            from = from, to = to,
                            group = group, `-log10(sig)` = pval, interaction_score = interaction_score
                        ),
                        aes(colour = group, alpha = `-log10(sig)`, width = interaction_score),
                        tension = 0.5
                    ) # + scale_edge_width(range = c(1, 3)) + scale_edge_alpha(limits = c(0, 1)) +
            } else {
                pl <- ggraph(gr, layout = "dendrogram", circular = TRUE) +
                    geom_conn_bundle(
                        data = get_con(
                            from = from, to = to,
                            group = group, `-log10(sig)` = pval, interaction_score = interaction_score
                        ),
                        aes(colour = group, alpha = interaction_score, width = `-log10(sig)`),
                        tension = 0.5
                    ) # + scale_edge_width(range = c(1, 3)) + scale_edge_alpha(limits = c(0, 1)) +
            }
            pl <- pl + scale_edge_color_manual(values = edge_group_colors) +
                geom_node_point(pch = 19, aes(
                    size = fraction, filter = leaf,
                    color = celltype, alpha = type
                )) + theme_void() + coord_fixed() +
                scale_size_continuous(limits = c(0, 1)) + scale_shape_manual(values = c(
                    ligand = 19,
                    receptor = 15
                )) + scale_color_manual(values = node_group_colors) +
                geom_text_repel(aes(x = x, y = y, label = label),
                    segment.square = TRUE,
                    segment.inflect = TRUE, segment.size = 0.2, force = 0.5,
                    size = 2, force_pull = 0
                ) + scale_alpha_manual(values = c(
                    ligand = 0.5,
                    receptor = 1
                )) + small_legend(keysize = 0.5) + ggtitle(input_group)
        } else {
            if (plot_score_as_thickness) {
                pl <- ggraph(gr, layout = "dendrogram", circular = TRUE) +
                    geom_conn_bundle(
                        data = get_con(
                            from = from, to = to,
                            `-log10(sig)` = pval, interaction_score = interaction_score
                        ),
                        aes(alpha = `-log10(sig)`, width = interaction_score),
                        tension = 0.5
                    )
            } else {
                pl <- ggraph(gr, layout = "dendrogram", circular = TRUE) +
                    geom_conn_bundle(
                        data = get_con(
                            from = from, to = to,
                            `-log10(sig)` = pval, interaction_score = interaction_score
                        ),
                        aes(alpha = interaction_score, width = `-log10(sig)`),
                        tension = 0.5
                    )
            }
            # scale_edge_width(range = c(1, 3)) +
            # scale_edge_alpha(limits = c(0, 1)) +
            pl <- pl + scale_edge_color_manual(values = edge_group_colors) +
                geom_node_point(pch = 19, aes(
                    size = fraction, filter = leaf,
                    color = celltype, alpha = type
                )) + theme_void() + coord_fixed() +
                scale_size_continuous(limits = c(0, 1)) + scale_shape_manual(values = c(
                    ligand = 19,
                    receptor = 15
                )) + scale_color_manual(values = node_group_colors) +
                geom_text_repel(aes(x = x, y = y, label = label),
                    segment.square = TRUE,
                    segment.inflect = TRUE, segment.size = 0.2, force = 0.5,
                    size = 2, force_pull = 0
                ) + # geom_node_text(aes(x = x*1.15, y=y*1.15, filter = leaf, label=label, size # =0.01)) + size
                scale_alpha_manual(values = c(ligand = 0.5, receptor = 1)) +
                small_legend(keysize = 0.5) + ggtitle(input_group)
        }
        return(pl)
    } else {
        return(NA)
    }
}

.chord_diagram4 <- function(tmp_dfx, lr_interactions, scaled, sep,
                            alpha, directional, show_legend, edge_cols, grid_cols, legend.pos.x, legend.pos.y,
                            title, grid_scale) {
    tmp_dfx <- .swap_ligand_receptor(tmp_dfx)
    unique_celltype <- unique(c(lr_interactions$`1`, lr_interactions$`2`))
    na_df <- data.frame(t(combn(unique_celltype, 2)))
    colnames(na_df) <- c("producer_swap", "receiver_swap")
    if (scaled) {
        interactions_items <- lr_interactions$scaled_means
    } else {
        interactions_items <- lr_interactions$means
    }
    names(interactions_items) <- paste0(lr_interactions$Var2, sep, lr_interactions$Var1)
    pvals_items <- lr_interactions$pvals
    names(pvals_items) <- paste0(lr_interactions$Var2, sep, lr_interactions$Var1)
    interactions_items[is.na(pvals_items)] <- 1
    tmp_dfx$pair_swap <- gsub("_", " - ", tmp_dfx$pair_swap)
    tmp_dfx$value <- interactions_items[tmp_dfx$barcode]
    tmp_dfx$pval <- pvals_items[tmp_dfx$barcode]
    edge_color <- .scPalette(length(unique(tmp_dfx$pair_swap)))
    names(edge_color) <- unique(tmp_dfx$pair_swap)
    if (!is.null(edge_cols)) {
        edge_color[names(edge_cols)] <- edge_cols
    }
    if (!is.null(grid_cols)) {
        if (length(grid_cols) != length(unique(tmp_dfx$receiver_swap))) {
            stop(paste0(
                "Please provide ", length(unique(tmp_dfx$receiver_swap)),
                " to grid_colors."
            ))
        } else {
            grid_color <- grid_cols
        }
    } else {
        grid_color <- .scPalette(length(unique(tmp_dfx$receiver_swap)))
    }
    if (is.null(grid_cols)) {
        names(grid_color) <- unique(tmp_dfx$receiver_swap)
    }
    tmp_dfx$edge_color <- edge_color[tmp_dfx$pair_swap]
    requireNamespace("colorspace")
    tmp_dfx$edge_color <- colorspace::adjust_transparency(tmp_dfx$edge_color,
        alpha = alpha
    )
    tmp_dfx$edge_color[is.na(tmp_dfx$pval)] <- NA
    tmp_dfx$grid_color <- grid_color[tmp_dfx$receiver_swap]
    tmp_dfx$grid_color[is.na(tmp_dfx$pval)] <- NA
    tmp_dfx <- tmp_dfx[!duplicated(tmp_dfx$barcode), ]
    # filter to non na
    tmp_dfx_not_na <- tmp_dfx[!is.na(tmp_dfx$pval), ]
    emptydf <- data.frame(matrix(ncol = ncol(tmp_dfx_not_na), nrow = nrow(na_df)))
    colnames(emptydf) <- colnames(tmp_dfx_not_na)
    emptydf$producer_swap <- na_df$producer_swap
    emptydf$receiver_swap <- na_df$receiver_swap
    tmp_dfx <- rbind(tmp_dfx_not_na, emptydf)
    tmp_dfx$value[is.na(tmp_dfx$value)] <- grid_scale
    if (directional == 2) {
        link.arr.type <- "triangle"
    } else {
        link.arr.type <- "big.arrow"
    }
    cells <- unique(c(tmp_dfx$producer_swap, tmp_dfx$receiver_swap))
    names(cells) <- cells
    circos.clear()
    chordDiagram(tmp_dfx[c("producer_swap", "receiver_swap", "value")],
        directional = directional,
        direction.type = c("diffHeight", "arrows"), link.arr.type = link.arr.type,
        annotationTrack = c("name", "grid"), col = tmp_dfx$edge_color, grid.col = grid_color,
        group = cells
    )
    requireNamespace("grid")
    requireNamespace("ComplexHeatmap")
    if (show_legend) {
        lgd <- ComplexHeatmap::Legend(
            at = names(edge_color), type = "grid",
            legend_gp = grid::gpar(fill = edge_color), title = "interactions"
        )
        ComplexHeatmap::draw(lgd,
            x = grid::unit(1, "npc") - grid::unit(legend.pos.x, "mm"),
            y = grid::unit(legend.pos.y, "mm"), just = c("right", "bottom")
        )
    }
    requireNamespace("graphics")
    graphics::title(main = title)
    circos.clear()
    gg <- recordPlot()
    return(gg)
}

.chord_diagram3 <- function(tmp_dfx, lr_interactions, scaled, sep,
                            alpha, directional, show_legend, edge_cols, grid_cols, legend.pos.x, legend.pos.y,
                            title) {
    tmp_dfx <- .swap_ligand_receptor(tmp_dfx)
    if (scaled) {
        interactions_items <- lr_interactions$scaled_means
    } else {
        interactions_items <- lr_interactions$means
    }
    names(interactions_items) <- paste0(lr_interactions$Var2, sep, lr_interactions$Var1)
    pvals_items <- lr_interactions$pvals
    names(pvals_items) <- paste0(lr_interactions$Var2, sep, lr_interactions$Var1)
    interactions_items[is.na(pvals_items)] <- 1
    tmp_dfx$pair_swap <- gsub("_", " - ", tmp_dfx$pair_swap)
    tmp_dfx$value <- interactions_items[tmp_dfx$barcode]
    tmp_dfx$pval <- pvals_items[tmp_dfx$barcode]
    edge_color <- .scPalette(length(unique(tmp_dfx$pair_swap)))
    names(edge_color) <- unique(tmp_dfx$pair_swap)
    if (!is.null(edge_cols)) {
        edge_color[names(edge_cols)] <- edge_cols
    }
    if (!is.null(grid_cols)) {
        if (length(grid_cols) != length(unique(tmp_dfx$receiver_swap))) {
            stop(paste0(
                "Please provide ", length(unique(tmp_dfx$receiver_swap)),
                " to grid_colors."
            ))
        } else {
            grid_color <- grid_cols
        }
    } else {
        grid_color <- .scPalette(length(unique(tmp_dfx$receiver_swap)))
    }
    if (is.null(grid_cols)) {
        names(grid_color) <- unique(tmp_dfx$receiver_swap)
    }
    tmp_dfx$edge_color <- edge_color[tmp_dfx$pair_swap]
    requireNamespace("colorspace")
    tmp_dfx$edge_color <- colorspace::adjust_transparency(tmp_dfx$edge_color,
        alpha = alpha
    )
    tmp_dfx$edge_color[is.na(tmp_dfx$pval)] <- NA
    tmp_dfx$grid_color <- grid_color[tmp_dfx$receiver_swap]
    tmp_dfx$grid_color[is.na(tmp_dfx$pval)] <- NA
    tmp_dfx <- tmp_dfx[!duplicated(tmp_dfx$barcode), ]
    if (directional == 2) {
        link.arr.type <- "triangle"
    } else {
        link.arr.type <- "big.arrow"
    }
    cells <- unique(c(tmp_dfx$producer_swap, tmp_dfx$receiver_swap))
    names(cells) <- cells
    circos.clear()
    chordDiagram(tmp_dfx[c("producer_swap", "receiver_swap", "value")],
        directional = directional,
        direction.type = c("diffHeight", "arrows"), link.arr.type = link.arr.type,
        annotationTrack = c("name", "grid"), col = tmp_dfx$edge_color, grid.col = grid_color,
        group = cells
    )
    requireNamespace("grid")
    requireNamespace("ComplexHeatmap")
    if (show_legend) {
        lgd <- ComplexHeatmap::Legend(
            at = names(edge_color), type = "grid",
            legend_gp = grid::gpar(fill = edge_color), title = "interactions"
        )
        ComplexHeatmap::draw(lgd,
            x = grid::unit(1, "npc") - grid::unit(legend.pos.x, "mm"),
            y = grid::unit(legend.pos.y, "mm"), just = c("right", "bottom")
        )
    }
    requireNamespace("graphics")
    graphics::title(main = title)
    circos.clear()
    gg <- recordPlot()
    return(gg)
}
