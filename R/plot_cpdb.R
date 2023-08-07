#' Plotting cellphonedb results
#'
#' @param cell_type1 cell type 1
#' @param cell_type2 cell type 2
#' @param scdata single-cell data. can be seurat/summarizedexperiment object
#' @param idents vector holding the idents for each cell or column name of scdata's metadata. MUST match cpdb's columns
#' @param means object holding means.txt from cpdb output
#' @param pvals object holding pvals.txt from cpdb output. Use relevant_interactions.txt if degs_analysis mode.
#' @param max_size max size of points.
#' @param keep_significant_only logical. Default is FALSE. Switch to TRUE if you only want to plot the significant hits from cpdb.
#' @param split.by column name in the metadata/coldata table to split the spots by. Can only take columns with binary options. If specified, name to split by MUST be specified in the meta file provided to cpdb prior to analysis.
#' @param gene.family default = NULL. some predefined group of genes. can take one (or several) of these default options: 'chemokines', 'Th1', 'Th2', 'Th17', 'Treg', 'costimulatory', 'coinhibitory', 'niche'. Also accepts name(s) of custom gene families.
#' @param custom_gene_family default = NULL. If provided, will update the gene.family function with this custom entry. Both `gene.family` (name of the custom family) and `custom_gene_family` must be specified for this to work. Provide either a data.frame with column names as name of family and genes in rows or a named likes like : list("customfamily" = c("genea", "geneb", "genec"))
#' @param genes default = NULL. can specify custom list of genes if gene.family is NULL
#' @param standard_scale logical. scale the expression to range from 0 to 1. NULL defaults to FALSE.
#' @param cluster_rows logical. whether or not to cluster the rows.
#' @param col_option specify plotting colours
#' @param default_style default = TRUE. Show all mean values and trace significant interactions with `higlight` colour. If FALSE, significant interactions will be presented as a white ring.
#' @param highlight_col colour for highlighting p <0.05
#' @param highlight_size stroke size for highlight if p < 0.05. if NULL, scales to -log10(pval).
#' @param special_character_regex_pattern search pattern if the cell type names contains special character. NULL defaults to '/|:|\\?|\\*|\\+|[\\]|\\(|\\)'.
#' @param degs_analysis if is cellphonedb degs_analysis mode.
#' @param return_table whether or not to return as a table rather than to plot.
#' @param exclude_interactions if provided, the interactions will be removed from the output.
#' @param ... passes arguments to grep for cell_type1 and cell_type2.
#' @return ggplot dot plot object of cellphone db output
#' @examples
#' \donttest{
#' data(kidneyimmune)
#' data(cpdb_output)
#' plot_cpdb("B cell", "CD4T cell", kidneyimmune, "celltype", means, pvals, split.by = "Experiment", genes = c("CXCL13", "CD274", "CXCR5"))
#' plot_cpdb("B cell", "CD4T cell", kidneyimmune, "celltype", means, pvals, split.by = "Experiment", gene.family = "chemokines")
#' }
#' @include utils.R
#' @import viridis
#' @import ggplot2
#' @import reshape2
#' @export

plot_cpdb <- function(
    cell_type1,
    cell_type2,
    scdata,
    idents,
    means,
    pvals,
    max_size = 8,
    keep_significant_only = FALSE,
    split.by = NULL,
    gene.family = NULL,
    custom_gene_family = NULL,
    genes = NULL,
    standard_scale = TRUE,
    cluster_rows = TRUE,
    col_option = viridis::viridis(50),
    default_style = TRUE,
    highlight_col = "red",
    highlight_size = NULL,
    special_character_regex_pattern = NULL,
    degs_analysis = FALSE,
    return_table = FALSE,
    exclude_interactions = NULL,
    ...) {
    requireNamespace("SingleCellExperiment")
    requireNamespace("grDevices")
    if (is.null(special_character_regex_pattern)) {
        special_character_regex_pattern <- DEFAULT_SPEC_PAT
    }
    if (class(scdata) %in% c("SingleCellExperiment", "SummarizedExperiment")) {
        metadata <- SingleCellExperiment::colData(scdata)
    } else if (class(scdata) == "Seurat") {
        metadata <- scdata@meta.data
    }
    means_mat <- .prep_table(means)
    pvals_mat <- .prep_table(pvals)
    if (degs_analysis) {
        col_start <- ifelse(colnames(pvals_mat)[13] == "classification", 14, 12)
        pvals_mat[, col_start:ncol(pvals_mat)] <- 1 - pvals_mat[, col_start:ncol(pvals_mat)]
    }
    cell_type1 <- .sub_pattern(cell_type = cell_type1, pattern = special_character_regex_pattern)
    cell_type2 <- .sub_pattern(cell_type = cell_type2, pattern = special_character_regex_pattern)
    if (is.null(genes)) {
        if (!is.null(gene.family)) {
            query_group <- .prep_query_group(data = means_mat, custom_gene_family = custom_gene_family)
            if (any(!gene.family %in% names(query_group))) {
                stop(paste("gene_family needs to be one of the following:", paste(names(query_group), collapse = ", ")))
            }
            query <- unique(do.call(c, query_group[tolower(gene.family)]))
        } else {
            query <- grep("", means_mat$interacting_pair, value = TRUE)
        }
    } else {
        if (!is.null(gene.family)) {
            stop("Please specify either genes or gene.family, not both.")
        } else {
            query <- grep(paste(genes, collapse = "|"), means_mat$interacting_pair, value = TRUE)
        }
    }
    # prepare the cell_type query
    if (!is.null(split.by)) {
        labels <- paste0(metadata[[split.by]], "_", metadata[[idents]])
        if (is.factor(metadata[[split.by]]) & is.factor(metadata[[idents]])) {
            labels <- factor(labels, levels = paste0(levels(metadata[[split.by]]), "_", rep(levels(metadata[[idents]]), each = length(levels(metadata[[split.by]])))))
        } else {
            labels <- factor(labels)
        }
        labels <- levels(labels)
        groups <- factor(metadata[[split.by]])
        groups <- levels(groups)
        if (length(groups) > 0) {
            # the purpose for this step is to allow for special characters to
            # be used in the celltype grepping
            if (length(groups) > 1) {
                labels2 <- gsub(paste0(paste0(groups, "_"), collapse = "|"), "", labels)
            } else {
                labels2 <- gsub(paste0(groups, "_"), "", labels)
            }
            # this returns the indices from the labels
            ct1 <- grep(cell_type1, labels2, ...)
            ct2 <- grep(cell_type2, labels2, ...)
            c_type1 <- as.list(labels[ct1])
            c_type2 <- as.list(labels[ct2])
            c_type1 <- lapply(c_type1, .sub_pattern, pattern = special_character_regex_pattern)
            c_type2 <- lapply(c_type2, .sub_pattern, pattern = special_character_regex_pattern)
            grp <- as.list(groups)
            celltype <- list()
            for (i in 1:length(c_type1)) {
                celltype[[i]] <- .create_celltype_query(ctype1 = c_type1[[i]], ctype2 = c_type2, sep = DEFAULT_SEP)
                celltype[[i]] <- lapply(grp, .keep_interested_groups, ct = celltype[[i]], sep = DEFAULT_SEP)
            }
            for (i in 1:length(celltype)) {
                celltype[[i]] <- celltype[[i]][-which(celltype[[i]] == "")]
            }
            celltype <- lapply(celltype, unlist)
            if (any(unlist(lapply(celltype, is.null)))) {
                rm <- which(unlist(lapply(celltype, is.null)))
                celltype <- celltype[-rm]
            }
            cell_type <- do.call(paste0, list(celltype, collapse = "|"))
        } else {
            labels <- metadata[[idents]]
            labels <- factor(labels)
            labels <- levels(labels)
            c_type1 <- as.list(grep(cell_type1, labels, value = TRUE, ...))
            c_type2 <- as.list(grep(cell_type2, labels, value = TRUE, ...))
            c_type1 <- lapply(c_type1, .sub_pattern, pattern = special_character_regex_pattern)
            c_type2 <- lapply(c_type2, .sub_pattern, pattern = special_character_regex_pattern)
            celltype <- list()
            for (i in 1:length(c_type1)) {
                celltype[[i]] <- .create_celltype_query(ctype1 = c_type1[[i]], ctype2 = c_type2, sep = DEFAULT_SEP)
            }
            cell_type <- do.call(paste0, list(celltype, collapse = "|"))
        }
    } else {
        labels <- metadata[[idents]]
        labels <- factor(labels)
        labels <- levels(labels)
        c_type1 <- as.list(grep(cell_type1, labels, value = TRUE))
        c_type2 <- as.list(grep(cell_type2, labels, value = TRUE))
        c_type1 <- lapply(c_type1, .sub_pattern, pattern = special_character_regex_pattern)
        c_type2 <- lapply(c_type2, .sub_pattern, pattern = special_character_regex_pattern)
        celltype <- list()
        for (i in 1:length(c_type1)) {
            celltype[[i]] <- .create_celltype_query(ctype1 = c_type1[[i]], ctype2 = c_type2, sep = DEFAULT_SEP)
        }
        cell_type <- do.call(paste0, list(celltype, collapse = "|"))
    }
    if (!is.null(gene.family) & is.null(genes)) {
        if (length(gene.family) == 1) {
            means_mat <- .prep_data_querygroup_celltype1(.data = means_mat, .query_group = query_group, .gene.family = gene.family, .cell_type = cell_type, .celltype = celltype, ...)
            pvals_mat <- .prep_data_querygroup_celltype1(.data = pvals_mat, .query_group = query_group, .gene.family = gene.family, .cell_type = cell_type, .celltype = celltype, ...)
        } else if (length(gene.family) > 1) {
            means_mat <- .prep_data_querygroup_celltype2(.data = means_mat, .query_group = query_group, .gene.family = gene.family, .cell_type = cell_type, .celltype = celltype, ...)
            pvals_mat <- .prep_data_querygroup_celltype2(.data = pvals_mat, .query_group = query_group, .gene.family = gene.family, .cell_type = cell_type, .celltype = celltype, ...)
        }
    } else if (is.null(gene.family) & !is.null(genes) | is.null(gene.family) & is.null(genes)) {
        means_mat <- .prep_data_query_celltype(.data = means_mat, .query = query, .cell_type = cell_type, .celltype = celltype, ...)
        pvals_mat <- .prep_data_query_celltype(.data = pvals_mat, .query = query, .cell_type = cell_type, .celltype = celltype, ...)
    }
    if (length(means_mat) == 0) {
        stop("Please check your options for split.by and your celltypes.")
    }
    # rearrange the columns so that it interleaves the two groups
    if (!is.null(split.by)) {
        if (length(groups) > 0) {
            grp <- as.list(groups)
            group_i <- lapply(grp, function(g) {
                gx <- grep(g, colnames(means_mat), ...)
                return(gx)
            })
            group_id <- do.call(c, group_i)
            means_mat <- means_mat[, as.vector(group_id), drop = FALSE]
            if (dim(pvals_mat)[2] > 0) {
                pvals_mat <- pvals_mat[, as.vector(group_id), drop = FALSE]
            } else {
                stop("No significant hits.")
            }
        }
    }
    if (keep_significant_only) {
        if (dim(pvals_mat)[2] == 0) {
            stop("No significant hits.")
        }
    }
    if (cluster_rows) {
        if (nrow(means_mat) > 2) {
            d <- dist(as.data.frame(means_mat))
            h <- hclust(d)
            means_mat <- means_mat[h$order, , drop = FALSE]
            pvals_mat <- pvals_mat[h$order, , drop = FALSE]
        }
    }
    # scaling
    if (standard_scale) {
        means_mat2 <- apply(means_mat, 1, range01)
        means_mat2 <- t(means_mat2)
    } else {
        means_mat2 <- means_mat
    }
    pvals_mat2 <- as.matrix(pvals_mat)
    means_mat2 <- as.matrix(means_mat2)
    xx <- which(means_mat == 0)
    if (length(xx) > 0) {
        means_mat2[which(means_mat == 0)] <- NA
    }
    # remove rows that are entirely NA
    pvals_mat2 <- pvals_mat2[rowSums(is.na(means_mat2)) != ncol(means_mat2), , drop = FALSE]
    means_mat2 <- means_mat2[rowSums(is.na(means_mat2)) != ncol(means_mat2), , drop = FALSE]
    if (standard_scale) {
        df_means <- melt(means_mat2, value.name = "scaled_means")
    } else {
        df_means <- melt(means_mat2, value.name = "means")
    }
    df_pvals <- melt(pvals_mat2, value.name = "pvals")
    df <- data.frame(cbind(df_means, pvals = df_pvals$pvals))
    xp <- which(df$pvals == 1)
    if (length(xp) > 0) {
        df$pvals[which(df$pvals == 1)] <- NA
    }
    if (keep_significant_only) {
        # keep the entire row/ all the comparisons
        df_ <- split(df, as.character(df$Var1))
        anysig <- lapply(df_, function(x) {
            keep <- any(x$pvals < 0.05)
            return(keep)
        })
        df_ <- df_[which(unlist(anysig))]
        names(df_) <- NULL
        df <- do.call(rbind, df_)
    }
    df$pvals[which(df$pvals == 0)] <- 0.001
    df$pvals[which(df$pvals >= 0.05)] <- NA

    if (!is.null(split.by)) {
        if (length(groups) > 0) {
            grp <- as.list(groups)
            grp2 <- lapply(grp, function(i) {
                x <- paste0(i, "_")
                return(x)
            })
            searchterm <- do.call(paste, list(grp2, collapse = "|"))
            df$group <- gsub(searchterm, "", df$Var2)
        }
    } else {
        df$group <- df$Var2
    }
    if (keep_significant_only) {
        if (standard_scale) {
            if (length(df$scaled_means) == 0) {
                stop("No significant genes found and plotting will not proceed.")
            }
        } else {
            if (length(df$means) == 0) {
                stop("No significant genes found and plotting will not proceed.")
            }
        }
    }
    df$Var2 <- gsub(DEFAULT_SEP, "-", df$Var2)
    final_levels <- unique(df$Var2)
    df$Var2 <- factor(df$Var2, unique(df$Var2))
    df$x_means_ <- df[, colnames(df_means)[3]]
    df$x_means_[df[, colnames(df)[4]] < 0.05] <- NA
    df$x_stroke <- df$x_means_
    df$x_stroke[!is.na(df$x_stroke)] <- 0
    df$x_stroke[is.na(df$x_stroke)] <- 2
    if (!is.null(exclude_interactions)) {
        df <- df[!df$Var1 %in% c(exclude_interactions), ]
    }
    df$significant <- ifelse(df$pvals < 0.05, "yes", NA)
    if (all(is.na(df$significant))) {
        df$significant <- "no"
        highlight_col <- "#ffffff"
    }
    if (return_table) {
        return(df)
    } else {
        if (default_style) {
            if (standard_scale) {
                g <- ggplot(df, aes(x = Var2, y = Var1, color = significant, fill = scaled_means, size = scaled_means))
            } else {
                g <- ggplot(df, aes(x = Var2, y = Var1, color = significant, fill = means, size = means))
            }
        } else {
            if (all(df$significant == "no")) {
                if (standard_scale) {
                    g <- ggplot(df, aes(x = Var2, y = Var1, color = significant, fill = scaled_means, size = scaled_means))
                } else {
                    g <- ggplot(df, aes(x = Var2, y = Var1, color = significant, fill = means, size = means))
                }
                default_style <- TRUE
            } else {
                highlight_col <- "#FFFFFF" # enforce this
                if (standard_scale) {
                    if (!is.null(highlight_size)) {
                        g <- ggplot(df, aes(x = Var2, y = Var1, fill = significant, colour = scaled_means, size = scaled_means, stroke = highlight_size))
                    } else {
                        g <- ggplot(df, aes(x = Var2, y = Var1, fill = significant, colour = scaled_means, size = scaled_means, stroke = x_stroke))
                    }
                } else {
                    if (!is.null(highlight_size)) {
                        g <- ggplot(df, aes(x = Var2, y = Var1, fill = significant, colour = means, size = means, stroke = highlight_size))
                    } else {
                        g <- ggplot(df, aes(x = Var2, y = Var1, fill = significant, colour = means, size = means, stroke = x_stroke))
                    }
                }
            }
        }
        g <- g + geom_point(pch = 21, na.rm = TRUE) +
            theme_bw() +
            theme(axis.text.x = element_text(angle = 45, hjust = 0, color = "#000000"), axis.text.y = element_text(color = "#000000")) +
            scale_x_discrete(position = "top") +
            scale_color_gradientn(colors = highlight_col, na.value = "white") +
            scale_size_continuous(range = c(0, max_size), aesthetics = "size") +
            scale_size_continuous(range = c(0, max_highlight_size), aesthetics = "stroke")
        if (default_style) {
            g <- g + scale_colour_manual(values = highlight_col, na.translate = FALSE) +
                guides(fill = guide_colourbar(barwidth = 4, label = TRUE, ticks = TRUE, draw.ulim = TRUE, draw.llim = TRUE, order = 1), size = guide_legend(reverse = TRUE, order = 2), stroke = guide_legend(reverse = TRUE, order = 3))
            if (length(col_option) == 1) {
                g <- g + scale_fill_gradientn(colors = grDevices::colorRampPalette(c("white", col_option))(100), na.value = "white")
            } else {
                g <- g + scale_fill_gradientn(colors = c("white", grDevices::colorRampPalette(col_option)(99)), na.value = "white")
            }
        } else {
            g <- g + scale_fill_manual(values = highlight_col, na.translate = FALSE) +
                guides(colour = guide_colourbar(barwidth = 4, label = TRUE, ticks = TRUE, draw.ulim = TRUE, draw.llim = TRUE, order = 1), size = guide_legend(reverse = TRUE, order = 2), stroke = guide_legend(reverse = TRUE, order = 3))

            df2 <- df
            if (standard_scale) {
                df2$scaled_means[df$pvals < 0.05] <- NA
                g <- g + geom_point(aes(x = Var2, y = Var1, colour = scaled_means, size = scaled_means), data = df2, inherit_aes = FALSE, na_rm = TRUE)
            } else {
                df2$means[df$pvals < 0.05] <- NA
                g <- g + geom_point(aes(x = Var2, y = Var1, colour = means, size = means), data = df2, inherit_aes = FALSE, na_rm = TRUE)
            }
            if (length(col_option) == 1) {
                g <- g + scale_colour_gradientn(colors = grDevices::colorRampPalette(c("white", col_option))(100), na.value = "white")
            } else {
                g <- g + scale_colour_gradientn(colors = c("white", grDevices::colorRampPalette(col_option)(99)), na.value = "white")
            }
        }
        if (!is.null(highlight_size)) {
            g <- g + guides(stroke = "none")
        }
        if (title != "") {
            g <- g + ggtitle(title)
        } else if (!is.null(gene.family) & is.null(genes)) {
            if (length(gene.family) > 1) {
                gene.family <- paste(gene.family, collapse = ", ")
            }
            g <- g + ggtitle(gene.family)
        }
        return(g)
    }
}
