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
    hues = seq(15, 375, length = n + 1)
    hcl(h = hues, l = 65, c = 100)[1:n]
}

.create_celltype_query <- function(ctype1, ctype2, sep) {
    ct1 = list()
    ct2 = list()
    for (i in 1:length(ctype2)) {
        ct1[i] = paste0("^", ctype1, sep, ctype2[i], "$")
        ct2[i] = paste0("^", ctype2[i], sep, ctype1, "$")
    }
    ct_1 = do.call(paste0, list(ct1, collapse = "|"))
    ct_2 = do.call(paste0, list(ct2, collapse = "|"))
    ct = list(ct_1, ct_2)
    ct = do.call(paste0, list(ct, collapse = "|"))
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
            sce_[which(SingleCellExperiment::rowData(sce_)[, "index"] %in% genes),
                ]
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

            sce_[which(SingleCellExperiment::rowData(sce_)[, "index"] %in% genes),
                ]
        }
    })
    cm <- mean(Matrix::rowMeans(SingleCellExperiment::counts(scex) > 0))
    return(cm)
}

.generateDf <- function(ligand, sep, receptor, receptor_a, receptor_b, pair, converted_pair,
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
    if (!is.null(splitted)){
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
                  x <- .cellTypeExpr_complex(sce_altx[[producers[i]]],
                    ligand[j], gsm)
                  y <- .cellTypeFraction_complex(sce_altx[[producers[i]]],
                    ligand[j], gsm)
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
                  x <- .cellTypeExpr_complex(sce_altx[[receivers[i]]],
                    receptor[j], gsm)
                  y <- .cellTypeFraction_complex(sce_altx[[receivers[i]]],
                    receptor[j], gsm)
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
            ra <- receptor_a[j]
            rb <- receptor_b[j]
            pr <- pair[j]
            out = data.frame(c(lg, rcp, ra, rb, pr, px, rx, producer_expression[lg,
                px], producer_fraction[lg, px], receiver_expression[rcp, rx], receiver_fraction[rcp,
                rx]))
            test_df = c(test_df, out)
        }
    }
    df_ <- do.call(rbind, test_df)
    row.names(df_) <- 1:nrow(df_)
    colnames(df_) <- c("ligand", "receptor", "receptor_a", "receptor_b", "pair",
        "producer", "receiver", "producer_expression", "producer_fraction", "receiver_expression",
        "receiver_fraction")
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
.swap_ligand_receptor <- function(df) {
    is_r_a = as.logical(df$receptor_a)
    is_r_b = as.logical(df$receptor_b)
    lg = df$ligand
    rp = df$receptor
    from = df$from
    to = df$to
    prd = df$producer
    rec = df$receiver
    prd_exp = df$producer_expression
    prd_fra = df$producer_fraction
    rec_exp = df$receiver_expression
    rec_fra = df$receiver_fraction
    # create swaps
    lg_swap = c()
    rp_swap = c()
    from_swap = c()
    to_swap = c()
    prd_swap = c()
    rec_swap = c()
    prd_exp_swap = c()
    prd_fra_swap = c()
    rec_exp_swap = c()
    rec_fra_swap = c()
    for (i in seq_along(is_r_a)) {
        if (!is_r_a[i]) {
            if (is_r_b[i]) {
                lg_swap = c(lg_swap, lg[i])
                rp_swap = c(rp_swap, rp[i])
                from_swap = c(from_swap, from[i])
                to_swap = c(to_swap, to[i])
                prd_swap = c(prd_swap, prd[i])
                rec_swap = c(rec_swap, rec[i])
                prd_exp_swap = c(prd_exp_swap, prd_exp[i])
                prd_fra_swap = c(prd_fra_swap, prd_fra[i])
                rec_exp_swap = c(rec_exp_swap, rec_exp[i])
                rec_fra_swap = c(rec_fra_swap, rec_fra[i])
            } else {
                lg_swap = c(lg_swap, lg[i])
                rp_swap = c(rp_swap, rp[i])
                from_swap = c(from_swap, from[i])
                to_swap = c(to_swap, to[i])
                prd_swap = c(prd_swap, prd[i])
                rec_swap = c(rec_swap, rec[i])
                prd_exp_swap = c(prd_exp_swap, prd_exp[i])
                prd_fra_swap = c(prd_fra_swap, prd_fra[i])
                rec_exp_swap = c(rec_exp_swap, rec_exp[i])
                rec_fra_swap = c(rec_fra_swap, rec_fra[i])
            }
        } else if (is_r_a[i]) {
            if (is_r_b[i]) {
                lg_swap = c(lg_swap, lg[i])
                rp_swap = c(rp_swap, rp[i])
                from_swap = c(from_swap, from[i])
                to_swap = c(to_swap, to[i])
                prd_swap = c(prd_swap, prd[i])
                rec_swap = c(rec_swap, rec[i])
                prd_exp_swap = c(prd_exp_swap, prd_exp[i])
                prd_fra_swap = c(prd_fra_swap, prd_fra[i])
                rec_exp_swap = c(rec_exp_swap, rec_exp[i])
                rec_fra_swap = c(rec_fra_swap, rec_fra[i])
            } else {
                lg_swap = c(lg_swap, rp[i])
                rp_swap = c(rp_swap, lg[i])
                from_swap = c(from_swap, to[i])
                to_swap = c(to_swap, from[i])
                prd_swap = c(prd_swap, rec[i])
                rec_swap = c(rec_swap, prd[i])
                prd_exp_swap = c(prd_exp_swap, rec_exp[i])
                prd_fra_swap = c(prd_fra_swap, rec_fra[i])
                rec_exp_swap = c(rec_exp_swap, prd_exp[i])
                rec_fra_swap = c(rec_fra_swap, prd_fra[i])
            }
        }
    }
    df$ligand_swap = lg_swap
    df$receptor_swap = rp_swap
    df$pair_swap = paste0(lg_swap, " - ", rp_swap)
    df$producer_swap = prd_swap
    df$receiver_swap = rec_swap
    df$producer_expression_swap = prd_exp_swap
    df$producer_fraction_swap = prd_fra_swap
    df$receiever_expression_swap = rec_exp_swap
    df$receiever_fraction_swap = rec_fra_swap
    df$from_swap = from_swap
    df$to_swap = to_swap
    return(df)
}