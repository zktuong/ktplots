#' Plotting genes as dotplot
#'
#' @param scdata single-cell data. can be seurat/summarizedexperiment object
#' @param idents column name holding the idents for each cell
#' @param genes genes you want to plot
#' @param split.by column name in the metadata/coldata table to split the spots by. If not provided, it will plot via idents provided.
#' @param pct.threshold float. required to keep gene expressed by minimum percentage of cells
#' @param scale logical. scale the expression to mean +/- SD. NULL defaults to TRUE.
#' @param standard_scale logical. scale the expression to range from 0 to 1. NULL defaults to FALSE.
#' @param keepLevels logical. keep the original factor of the levels of the idents (for plotting)
#' @param save.plot logical. will try to save the pdf
#' @param h height of plot
#' @param w width of plot
#' @param filepath path to file, or path to folder
#' @param filename path to file
#' @param heat_cols colour gradient for the dot plot
#' @param col_limits set limits to the color gradient
#' @param outline_col colour of outlines if fill = TRUE
#' @param outline_size stroke size of outlines if fill = TRUE
#' @return ggplot dot plot object of selected genes
#' @examples
#' \donttest{
#' data(kidneyimmune)
#' geneDotPlot(kidneyimmune, genes = c("CD68", "CD80", "CD86", "CD74", "CD2", "CD5"), idents = "celltype", split.by = "Project", standard_scale = TRUE) + theme(strip.text.x = element_text(angle = 45, hjust = 0))
#' }
#' @import dplyr
#' @import gtools
#' @import Matrix
#' @import ggplot2
#' @import reshape2
#' @import RColorBrewer
#' @export
geneDotPlot <- function(scdata, idents, genes, split.by = NULL, pct.threshold = 0.05, scale = NULL, standard_scale = NULL, keepLevels = TRUE, save.plot = FALSE, h = 5, w = 5, filepath = NULL, filename = NULL, heat_cols = NULL, col_limits = NULL, fill = FALSE, outline_col = "black", outline_size = .2) {
    if (class(scdata) %in% c("SingleCellExperiment", "SummarizedExperiment")) {
        cat("data provided is a SingleCellExperiment/SummarizedExperiment object", sep = "\n")
        cat("extracting expression matrix", sep = "\n")
        requireNamespace("SummarizedExperiment")
        exp_mat <- SummarizedExperiment::assay(scdata)
        metadata <- SummarizedExperiment::colData(scdata)
    } else if (class(scdata) == "Seurat") {
        cat("data provided is a Seurat object", sep = "\n")
        cat("extracting expression matrix", sep = "\n")
        requireNamespace("Seurat")
        exp_mat <- tryCatch(scdata@data, error = function(e) {
            tryCatch(Seurat::GetAssayData(object = scdata), error = function(e) {
                stop(paste0("are you sure that your data is normalized?"))
                return(NULL)
            })
        })
        metadata <- scdata@meta.data
    }

    cat(paste0("attempting to subset the expression matrix to the ", length(genes), " genes provided"), sep = "\n")
    # expr_mat_filtered <- exp_mat[row.names(exp_mat) %in% genes, ]
    # exp_mat <- as.matrix(exp_mat)
    expr_mat_filtered <- exp_mat[match(rev(genes), row.names(exp_mat))[!is.na(match(rev(genes), row.names(exp_mat)))], , drop = FALSE]

    cat(paste0("found ", dim(expr_mat_filtered)[1], " genes in the expression matrix", sep = "\n"))

    if (!is.null(split.by)) {
        labels <- paste0(as.character(metadata[[split.by]]), "_", as.character(metadata[[idents]]))
        labels <- factor(labels)
    } else {
        cat("no groups information provided. defaulting to idents only", sep = "\n")
        labels <- factor(metadata[[idents]])
    }

    cat("preparing the final dataframe ...", sep = "\n")
    quick_prep <- function(expr, label, groups. = NULL, scale. = scale, meta = metadata, id = idents, standard_scale. = standard_scale) {
        expr.df <- tryCatch(data.frame(label = label, t(as.matrix(expr)), check.names = FALSE), error = function(e) {
            data.frame(label = label, t(Matrix::Matrix(expr, sparse = FALSE)), check.names = FALSE)
        })

        meanExpr <- split(expr.df, expr.df$label)
        meanExpr <- lapply(meanExpr, function(x) {
            x <- x[, -1, drop = FALSE]
            x <- x %>% colMeans()
            return(x)
        })

        # names(meanExpr) <- levels(label)
        meanExpr <- do.call(rbind, meanExpr)

        # control the scaling here
        if (length(standard_scale.) > 0) {
            if (standard_scale.) {
                meanExpr_ <- apply(meanExpr, 2, range01)
            } else {
                meanExpr_ <- meanExpr
            }
        }

        if (length(scale.) < 1) {
            if (length(standard_scale.) > 0) {
                if (standard_scale.) {
                    meanExpr <- meanExpr_
                } else {
                    meanExpr <- meanExpr
                }
            } else {
                meanExpr <- scale(meanExpr)
            }
        } else {
            if (scale.) {
                if (length(standard_scale.) > 0) {
                    if (standard_scale.) {
                        meanExpr <- meanExpr_
                    } else {
                        meanExpr <- scale(meanExpr)
                    }
                } else {
                    meanExpr <- scale(meanExpr)
                }
            } else {
                meanExpr <- meanExpr
            }
        }

        label.list <- as.list(levels(label))
        names(label.list) <- levels(label)
        exp <- lapply(label.list, function(x) {
            exp_f <- expr.df %>%
                dplyr::filter(label == x) %>%
                dplyr::select(-matches("label"))
            return(exp_f)
        })

        cellNumbers <- do.call(rbind, lapply(exp, dim))[, 1]

        pct <- list()
        pct <- lapply(exp, function(y) sapply(y, function(x) length(which(x > 0))))

        names(pct) <- levels(label)
        pct <- do.call(rbind, pct)
        final.pct <- pct / cellNumbers

        meltedMeanExpr <- reshape2::melt(meanExpr)
        meltedfinal.pct <- reshape2::melt(final.pct)

        if (!is.null(groups)) {
            meltedMeanExpr$Var3 <- gsub(".*_", "", meltedMeanExpr$Var1)
            meltedfinal.pct$Var3 <- gsub(".*_", "", meltedfinal.pct$Var1)
            meltedfinal.pct <- meltedfinal.pct[order(meltedfinal.pct$Var3, meltedfinal.pct$Var2), ]
            meltedMeanExpr <- meltedMeanExpr[order(meltedMeanExpr$Var3, meltedMeanExpr$Var2), ]
            meltedfinal.pct <- meltedfinal.pct[, -4]
            meltedMeanExpr <- meltedMeanExpr[, -4]
        } else {
            meltedfinal.pct <- meltedfinal.pct[order(meltedfinal.pct$Var1, meltedfinal.pct$Var2), ]
            meltedMeanExpr <- meltedMeanExpr[order(meltedMeanExpr$Var1, meltedMeanExpr$Var2), ]
        }

        df <- cbind(meltedMeanExpr, meltedfinal.pct$value)
        if ((length(scale.) > 0 && scale.) | (length(scale.) < 1 && length(standard_scale.) < 1) | (length(standard_scale.) > 0 && standard_scale.)) {
            colnames(df) <- c("celltype", "gene", "scale.mean", "pct")
        } else {
            colnames(df) <- c("celltype", "gene", "mean", "pct")
        }

        # add some groupings
        if (!is.null(groups.)) {
            df$group <- groups.[1]
            for (i in 2:length(groups.)) {
                df$group[grep(groups.[i], df$celltype)] <- groups.[i]
            }
            df$group <- factor(df$group, levels = groups.)
            remove.pattern <- paste0(groups., "_", collapse = "|")
            df$cell_type <- gsub(pattern = remove.pattern, "", df$celltype)
            df$cell_type <- factor(df$cell_type, levels = levels(meta[[id]]))
            df <- df[with(df, order(df$cell_type, df$group)), ]
        } else {
            df$cell_type <- as.factor(df$celltype)
            df$group <- as.factor(df$celltype)
            df <- df[order(df$cell_type), ]
        }
        return(df)
    }

    if (!is.null(split.by)) {
        plot.df <- quick_prep(expr_mat_filtered, labels, levels(droplevels(factor(metadata[[split.by]]))))
    } else {
        plot.df <- quick_prep(expr_mat_filtered, labels)
    }

    if (!is.null(pct.threshold)) {
        cat(paste0("setting minimum percentage of cells expressing gene to be ", pct.threshold * 100, "% of cluster/cell-type"), sep = "\n")

        filter <- split(plot.df, plot.df$gene)
        remove.genes <- lapply(filter, function(x) {
            if (max(x$pct) < pct.threshold) {
                return(as.character(x$gene))
            }
        })
        remove.genes <- unique(unlist(remove.genes))
        keep.genes <- lapply(filter, function(x) {
            if (max(x$pct) >= pct.threshold) {
                return(as.character(x$gene))
            }
        })
        keep.genes <- unique(unlist(keep.genes))
        cat("the following genes are removed", sep = "\n")
        print(remove.genes)
    } else if (is.null(pct.threshold) | pct.threshold == 0) {
        warning("are you sure you don't want to set a cut off?")
        keep.genes <- plot.df %>%
            dplyr::select(gene) %>%
            unique() %>%
            unlist() %>%
            as.character()
    }

    # finally keep the genes for plotting
    plot.df.final <- plot.df[plot.df$gene %in% keep.genes, ]
    if (!keepLevels) {
        plot.df.final$cell_type <- factor(plot.df.final$cell_type, levels = gtools::mixedsort(levels(plot.df.final$cell_type)))
    } else {
        plot.df.final$cell_type <- plot.df.final$cell_type
    }

    plot.df.final$pct[plot.df.final$pct == 0] <- NA

    if (!is.null(heat_cols)) {
        heat_cols <- heat_cols
    } else {
        heat_cols <- rev(brewer.pal(9, "RdBu"))
    }

    # subset the plotting objects
    doplot <- function(obj, group. = NULL, file_name = filename, file_path = filepath, dim_w = w, dim_h = h, limits. = col_limits, do.plot = save.plot, scale. = scale, standard_scale. = standard_scale) {
        if (is.null(group.)) {
            if ((length(scale.) > 0 && scale.) | (length(scale.) < 1 && length(standard_scale.) < 1) | (length(standard_scale.) > 0 && standard_scale.)) {
                g <- ggplot(obj, aes(x = 0, y = gene, size = pct, colour = scale.mean))
            } else {
                g <- ggplot(obj, aes(x = 0, y = gene, size = pct, colour = mean))
            }
            g <- g + geom_point(pch = 16, na.rm = TRUE) +
                scale_y_discrete(position = "left") +
                scale_x_discrete(position = "bottom") +
                scale_colour_gradientn(colors = heat_cols, limits = limits., na.value = "grey90", oob = scales::squish) +
                scale_radius(range = c(0, 4), limits = c(0, 1)) +
                theme_bw() +
                theme(
                    axis.text.x = element_text(angle = 90, hjust = 1),
                    axis.title.x = element_blank(),
                    axis.ticks = element_blank(),
                    axis.title.y = element_blank(),
                    axis.line = element_blank(),
                    panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),
                    panel.border = element_blank(),
                    strip.background = element_blank()
                ) +
                facet_grid(~cell_type)
        } else {
            if ((length(scale.) > 0 && scale.) | (length(scale.) < 1 && length(standard_scale.) < 1) | (length(standard_scale.) > 0 && standard_scale.)) {
                g <- ggplot(obj, aes(x = group, y = gene, size = pct, colour = scale.mean))
            } else {
                g <- ggplot(obj, aes(x = group, y = gene, size = pct, colour = mean))
            }

            g <- g + geom_point(pch = 16, na.rm = TRUE) +
                scale_y_discrete(position = "left") +
                scale_x_discrete(position = "bottom") +
                scale_colour_gradientn(colors = heat_cols, limits = limits., na.value = "grey90", oob = scales::squish) +
                scale_radius(range = c(0, 4), limits = c(0, 1)) +
                theme_bw() +
                theme(
                    axis.text.x = element_text(angle = 90, hjust = 1),
                    axis.title.x = element_blank(),
                    axis.ticks = element_blank(),
                    axis.title.y = element_blank(),
                    axis.line = element_blank(),
                    panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),
                    panel.border = element_blank(),
                    strip.background = element_blank()
                ) +
                facet_grid(~cell_type)
        }

        if (do.plot) {
            if (is.null(file_name) && is.null(file_path)) {
                out_path <- "./geneDotPlot.df"
                warning("no file name provided. saving plot to ", getwd(), "/geneDotPlot.pdf")
                ggsave("./geneDotPlot.pdf", plot = g, width = dim_w, height = dim_h, device = "pdf", useDingbats = FALSE)
            } else if (!is.null(file_name) && is.null(file_path)) {
                cat(paste0("saving plot to ", file_name), sep = "\n")
                tryCatch(ggsave(file_name, plot = g, width = dim_w, height = dim_h, device = "pdf", useDingbats = FALSE), error = function(e) {
                    ggsave("./geneDotPlot.df", plot = g, width = dim_w, height = dim_h, device = "pdf", useDingbats = FALSE)
                    warning("file name provided is not suitable. saving as geneDotPlot.pdf")
                })
            } else if (is.null(file_name) && !is.null(file_path)) {
                cat(paste0("saving plot to ", file_path), sep = "\n")
                if (grepl(".pdf", file_path)) {
                    ggsave(file_path, plot = g, width = dim_w, height = dim_h, device = "pdf", useDingbats = FALSE)
                } else {
                    dir.create(file_path, recursive = TRUE)
                    ggsave(paste0(file_path, "/geneDotPlot.df"), plot = g, width = dim_w, height = dim_h, device = "pdf", useDingbats = FALSE)
                    warning(paste0("file path provided is not suitable. saving as ", file_path, "/geneDotPlot.pdf"))
                }
            } else if (!is.null(file_name) && !is.null(file_path)) {
                cat(paste0("saving plot to ", paste0(file_path, "/", file_name)), sep = "\n")
                dir.create(file_path, recursive = TRUE)
                tryCatch(ggsave(paste0(file_path, "/", file_name), plot = g, width = dim_w, height = dim_h, device = "pdf", useDingbats = FALSE), error = function(e) {
                    ggsave("./geneDotPlot.df", plot = g, width = dim_w, height = dim_h, device = "pdf", useDingbats = FALSE)
                    warning("file path provided is not suitable. saving as geneDotPlot.pdf")
                })
            }
        }
        return(g)
    }

    fillplot <- function(obj, group. = NULL, file_name = filename, file_path = filepath, dim_w = w, dim_h = h, limits. = col_limits, do.plot = save.plot, scale. = scale, standard_scale. = standard_scale, outline_col. = outline_col, outline_size. = outline_size) {
        if (is.null(group.)) {
            if ((length(scale.) > 0 && scale.) | (length(scale.) < 1 && length(standard_scale.) < 1) | (length(standard_scale.) > 0 && standard_scale.)) {
                g <- ggplot(obj, aes(x = 0, y = gene, size = pct, fill = scale.mean))
            } else {
                g <- ggplot(obj, aes(x = 0, y = gene, size = pct, fill = mean))
            }
            g <- g + geom_point(pch = 21, color = outline_col., stroke = outline_size.) +
                scale_y_discrete(position = "left") +
                scale_x_discrete(position = "bottom") +
                scale_fill_gradientn(colors = heat_cols, limits = limits., na.value = "grey90", oob = scales::squish) +
                scale_radius(range = c(0, 4), limits = c(0, 1)) +
                theme_bw() + theme(
                    axis.text.x = element_text(angle = 90, hjust = 1),
                    axis.title.x = element_blank(),
                    axis.ticks = element_blank(),
                    axis.title.y = element_blank(),
                    axis.line = element_blank(),
                    panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),
                    panel.border = element_blank(),
                    strip.background = element_blank()
                ) +
                facet_grid(~cell_type)
        } else {
            if ((length(scale.) > 0 && scale.) | (length(scale.) <
                1 && length(standard_scale.) < 1) | (length(standard_scale.) >
                0 && standard_scale.)) {
                g <- ggplot(obj, aes(
                    x = group, y = gene, size = pct,
                    fill = scale.mean
                ))
            } else {
                g <- ggplot(obj, aes(
                    x = group, y = gene, size = pct,
                    fill = mean
                ))
            }
            g <- g + geom_point(pch = 21, color = outline_col., stroke = outline_size.) + scale_y_discrete(position = "left") +
                scale_x_discrete(position = "bottom") +
                scale_fill_gradientn(colors = heat_cols, limits = limits., na.value = "grey90", oob = scales::squish) +
                scale_radius(range = c(0, 4), limits = c(0, 1)) +
                theme_bw() + theme(
                    axis.text.x = element_text(angle = 90, hjust = 1),
                    axis.title.x = element_blank(),
                    axis.ticks = element_blank(),
                    axis.title.y = element_blank(),
                    axis.line = element_blank(),
                    panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),
                    panel.border = element_blank(),
                    strip.background = element_blank()
                ) +
                facet_grid(~cell_type)
        }
        if (do.plot) {
            if (is.null(file_name) && is.null(file_path)) {
                out_path <- "./geneDotPlot.df"
                warning("no file name provided. saving plot to ", getwd(), "/geneDotPlot.pdf")
                ggsave("./geneDotPlot.pdf", plot = g, width = dim_w, height = dim_h, device = "pdf", useDingbats = FALSE)
            } else if (!is.null(file_name) && is.null(file_path)) {
                cat(paste0("saving plot to ", file_name), sep = "\n")
                tryCatch(ggsave(file_name, plot = g, width = dim_w, height = dim_h, device = "pdf", useDingbats = FALSE), error = function(e) {
                    ggsave("./geneDotPlot.df", plot = g, width = dim_w, height = dim_h, device = "pdf", useDingbats = FALSE)
                    warning("file name provided is not suitable. saving as geneDotPlot.pdf")
                })
            } else if (is.null(file_name) && !is.null(file_path)) {
                cat(paste0("saving plot to ", file_path), sep = "\n")
                if (grepl(".pdf", file_path)) {
                    ggsave(file_path, plot = g, width = dim_w, height = dim_h, device = "pdf", useDingbats = FALSE)
                } else {
                    dir.create(file_path, recursive = TRUE)
                    ggsave(paste0(file_path, "/geneDotPlot.df"), plot = g, width = dim_w, height = dim_h, device = "pdf", useDingbats = FALSE)
                    warning(paste0("file path provided is not suitable. saving as ", file_path, "/geneDotPlot.pdf"))
                }
            } else if (!is.null(file_name) && !is.null(file_path)) {
                cat(paste0("saving plot to ", paste0(file_path, "/", file_name)), sep = "\n")
                dir.create(file_path, recursive = TRUE)
                tryCatch(ggsave(paste0(file_path, "/", file_name), plot = g, width = dim_w, height = dim_h, device = "pdf", useDingbats = FALSE), error = function(e) {
                    ggsave("./geneDotPlot.df", plot = g, width = dim_w, height = dim_h, device = "pdf", useDingbats = FALSE)
                    warning("file path provided is not suitable. saving as geneDotPlot.pdf")
                })
            }
        }
        return(g)
    }
    if (fill) {
        gg <- fillplot(plot.df.final, split.by, file_name = filename, file_path = filepath, dim_w = w, dim_h = h, limits. = col_limits, do.plot = save.plot, scale. = scale, standard_scale. = standard_scale, outline_col. = outline_col, outline_size. = outline_size)
    } else {
        gg <- doplot(plot.df.final, split.by, file_name = filename, file_path = filepath, dim_w = w, dim_h = h, limits. = col_limits, do.plot = save.plot, scale. = scale, standard_scale. = standard_scale)
    }

    gg
    return(gg)
}
