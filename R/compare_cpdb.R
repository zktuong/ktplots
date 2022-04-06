#' Plotting cellphonedb results
#'
#' @param cpdb_meta data.frame containing the sample name, cellphonedb 'out' folder path (containing means.txt and pvalues.txt/relevant_interactions.txt if version 3), and single-cell object file path (.h5ad or .rds)
#' @param sample_metadata data.frame containing the sample name, and groupings.
#' @param celltypes subset celltypes for comparison
#' @param celltype_col column name in single cell object holding celtype annotation.
#' @param groupby for significance testing. only if method is t.test or wilcoxon.
#' @param formula for signfiicance testing. only if method is lme.
#' @param method one of 't.test', 'wilcox', 'lme'
#' @param BPPARAM BiocParallelParam class.
#' @param version3 boolean. if cellphonedb version 3
#' @param verbose Whether or not to print messages.
#' @param p.adjust.mode whether or not to correct all interactions, or to perform it on a celltype-celltype basis.
#' @param p.adjust.method passed to `method`in `stats::p.adjust`.
#' @param ... passed to tests.
#' @return results for plotting
#' @examples
#' \donttest{}
#' @import BiocParallel
#' @import dplyr
#' @export
compare_cpdb <- function(cpdb_meta, sample_metadata, celltypes, celltype_col, groupby = NULL, formula = NULL, method = c('t.test', 'wilcox', 'lme'), BPPARAM = SerialParam(), version3 = FALSE, verbose = TRUE, p.adjust.mode = c('celltype', 'all'), p.adjust.method = 'fdr', ...) {
    options(warn = -1)
    sample <- cpdb_meta[, 1]
    cpdb_out_folder <- cpdb_meta[, 2]
    expression_file <- cpdb_meta[, 3]
    names(cpdb_out_folder) <- sample
    names(expression_file) <- sample
    if (length(method) == 3){
        method = 't.test'
    }
    if (length(p.adjust.mode) == 2){
        method = 'celltype'
    }
    if (verbose) {
        cat("Reading cellphonedb outputs", sep = "\n")
    }
    means <- bplapply(cpdb_out_folder, function(x) read.delim(paste0(x, "/means.txt"),
        check.names = FALSE), BPPARAM = SerialParam(progressbar = verbose))
    if (version3){
        pval_file <- '/relevant_interactions.txt'
    } else {
        pval_file <- '/pvalues.txt'
    }
    pvals <- bplapply(cpdb_out_folder, function(x) read.delim(paste0(x, pval_file),
        check.names = FALSE), BPPARAM = SerialParam(progressbar = verbose))

    if (verbose) {
        cat("Reading single cell object", sep = "\n")
    }

    requireNamespace('Matrix')
    sces <- bplapply(expression_file, function(x) {
        if (endsWith(x, ".h5ad")) {
            require(reticulate)
            ad = import("anndata")
            adata = ad$read_h5ad(paste0(x))
            counts <- Matrix::t(adata$X)
            row.names(counts) <- row.names(adata$var)
            colnames(counts) <- row.names(adata$obs)
            sce <- SingleCellExperiment(list(counts = counts), colData = adata$obs,
                rowData = adata$var)
        } else {
            sce <- readRDS(x)
            if (class(sce) == "Seurat") {
                sce <- as.SingleCellExperiment(sce)
            }
        }
        return(sce)
    }, BPPARAM = SerialParam(progressbar = verbose))

    # set up cell type combinations
    combs <- t(combn(celltypes, 2))
    # add in self-self interactions too
    combs <- rbind(combs, cbind(as.character(celltypes), as.character(celltypes)))
    comb_list <- apply(combs, 1, as.list)

    get_means <- function(x) {
        if (!is.na(x)) {
            out <- x$means
            names(out) <- paste0(x$group, ">@<", x$Var1)
            if (any(is.na(out))) {
                out <- out[!is.na(out)]
            }
            return(out)
        } else {
            return(NA)
        }
    }

    # extract interactions with plot_cpdb
    if (verbose) {
        cat(paste0("Extracting cpdb results across ", dim(combs)[1], " celltype pairwise combinations"),
            sep = "\n")
    }

    ct_sigs <- bplapply(comb_list, function(x){
        bpmapply(function(sce, mean, pval) {
                s <- tryCatch(plot_cpdb(cell_type1 = paste0(x[1], '$'),
                    cell_type2 = paste0(x[2], '$'),
                    scdata = sce,
                    idents = celltype_col, # column name where the cell ids are located in the metadata
                    means = mean,
                    pvals = pval,
                    keep_significant_only = TRUE,
                    scale = FALSE,
                    return_table = TRUE,
                ), error = function(e) return(NA))
                out <- suppressWarnings(tryCatch(get_means(s), error = function(e) return(NA)))
                return(out)
            }, sces, means, pvals, BPPARAM = BPPARAM, SIMPLIFY = FALSE)
    }, BPPARAM = SerialParam(progressbar = verbose))

    # collapse to matrix
    requireNamespace("plyr")
    meta <- bplapply(ct_sigs, function(y) {
        z <- as.data.frame(plyr::rbind.fill(lapply(y, function(y) {
            tryCatch(as.data.frame(t(y), stringsAsFactors = FALSE), error = function(e) return(NA))
        })))
        row.names(z) <- names(sces)
        return(z)
    }, BPPARAM = BPPARAM)

    if (verbose) {
        cat("Filtering interactions", sep = "\n")
    }
    # remove interactions that are all empty
    meta2 <- bplapply(meta, function(y) {
        if (!is.logical(y)) {
            f <- apply(y, 2, function(x) all(is.na(x)))
            y <- y[, !f]
        } else {
            y <- NA
        }
        return(y)
    }, BPPARAM = BPPARAM)
    meta2 <- bplapply(meta2, function(x) {
        x[is.na(x)] <- 0
        return(x)
    }, BPPARAM = SerialParam(progressbar = verbose))

    if (verbose) {
        cat("Preparing data", sep = "\n")
    }
    # split each column into a list
    res2 <- bplapply(meta2, function(x) {
        tmp <- as.list(x)
        tmp <- bplapply(tmp, function(z) {
            z <- as.numeric(z)
            names(z) <- row.names(x)
            return(z)
        }, BPPARAM = BPPARAM)
        return(tmp)
    }, BPPARAM = SerialParam(progressbar = verbose))

    requireNamespace('reshape2')
    test_fun <- function(int_score, data, col, method, ...) {
        # 'C'
        data$int_score <- int_score
        if (class(data[, col]) == "factor") {
            lvls <- levels(data[, col])
        } else {
            data[, col] <- factor(data[, col])
            lvls <- levels(data[, col])
        }

        if (method == "wilcox") {
            test <- pairwise.wilcox.test(data$int_score, data[, col], p.adjust.method = "none", ...)
        } else if (method == "t.test") {
            test <- pairwise.t.test(data$int_score, data[, col], p.adjust.method = "none", ...)
        }

        tmp <- reshape2::melt(test$p.value)
        p <- tmp$value
        names(p) <- paste0(tmp$Var1, ">:<_vs_>:<", tmp$Var2)
        p <- p[!is.na(p)]
        diff <- group_by(data, get(col)) %>% summarise(mean = log1p(mean(int_score,
            na.rm = TRUE)))
        diff <- as.data.frame(diff)
        row.names(diff) <- diff[, 1]
        diff <- diff[, "mean", drop = FALSE]
        n <- names(p)
        n <- strsplit(n, ">:<_vs_>:<")
        lfc <- lapply(n, function(x) {
            diff[x[1], ] - diff[x[2], ]
        })
        names(lfc) <- gsub(">:<_vs_>:<", "_vs_", names(p))
        names(p) <- gsub(">:<_vs_>:<", "_vs_", names(p))
        names(p) <- paste0("P_", names(p))
        outdf <- cbind(t(data.frame(lfc)), data.frame(p))
        colnames(outdf) <- c("LFC", "pval")
        outdf$contrast <- rownames(outdf)
        return(outdf)
    }

    if (method != 'lme'){
        if (is.null(groupby)) {
            stop("Please provide column name for contrasts to be extracted.")
        }
        if (verbose) {
            cat(paste0("Running pairwise ", method), sep = "\n")
        }
        res3 <- bplapply(res2, function(x) {
            tmp <- bplapply(x, function(y) test_fun(int_score = y, data = sample_metadata,
                col = groupby, method = method, ...), BPPARAM = BPPARAM)
            tmp <- plyr::ldply(tmp, data.frame)
            return(tmp)
        }, BPPARAM = SerialParam(progressbar = verbose))
        res3 <- do.call(rbind, res3)
        tmpct <- as.data.frame(do.call(rbind, strsplit(res3[,1], '>@<'))[,1:2])
        tmpct[,3] <- paste0(tmpct[,1], '>@<', tmpct[,2])
        res3$celltypes <- tmpct[,3]
        res3 <- split(res3, res3$contrast)
    } else if (method == 'lme'){
        require(lmerTest)
        if (!is.null(formula)) {
            fullFormula <- update.formula(formula, int_score ~ ., simplify = FALSE,
                )
        } else {
            stop("Please provide the formula")
        }

        if (verbose) {
            cat("Running lmer model", sep = "\n")
        }
        res3_ <- bplapply(res2, function(x) {
            bplapply(x, function(y) {
                sample_metadata[, "int_score"] <- y
                results <- suppressMessages(suppressWarnings(tryCatch(lmer(fullFormula,
                    data = sample_metadata, ...), error = function(e) return(NA))))
                return(results)
            }, BPPARAM = BPPARAM)
        }, BPPARAM = SerialParam(progressbar = verbose))

        if (verbose) {
            cat("Extracting statistics", sep = "\n")
        }

        res3fitstats <- bplapply(res3_, function(fit) {
            if (length(fit) > 0) {
                if (!is.na(fit)) {
                    out <- lapply(fit, function(fit2) {
                        if (!is.na(fit2)) {
                          y <- suppressWarnings(summary(fit2))
                          E <- y$coefficients[-1, 1]
                          if (length(E) > 1) {
                            names(E) <- paste0("fit_estimates_", names(E))
                          } else {
                            names(E) <- paste0("fit_estimates_", rownames(y$coefficients)[2])
                          }
                          P <- y$coefficients[-1, 5]
                          if (length(E) > 1) {
                            names(P) <- paste0("fit_P_", names(P))
                          } else {
                            names(P) <- paste0("fit_P_", rownames(y$coefficients)[2])
                          }
                          return(c(E, P))
                        } else {
                          return(NA)
                        }
                    })
                    Singular <- lapply(fit, function(fit2) {
                        if (!is.na(fit2)) {
                          out <- as.numeric(isSingular(fit2))
                          names(out) <- "Singular"
                          return(out)
                        } else {
                          return(NA)
                        }
                    })
                    Conv <- lapply(fit, function(fit2) {
                        if (!is.na(fit2)) {
                          out <- length(slot(fit2, "optinfo")$conv$lme4$messages)
                          names(out) <- "Conv"
                          return(out)
                        } else {
                          return(NA)
                        }
                    })
                    anv <- lapply(fit, function(fit2) {
                        if (!is.na(fit2)) {
                          anva <- anova(fit2, ddf="Satterthwaite")
                          wp <- anva[, 6]
                          names(wp) <- paste0("anova_P_", row.names(anva))
                          return(wp)
                        } else {
                          return(NA)
                        }
                    })
                    return(list(anv = anv, fit = out, Singular = Singular, Conv = Conv))
                } else {
                    return(NA)
                }
            } else {
                return(NA)
            }
        }, BPPARAM = SerialParam(progressbar = verbose))
        res3fitstats2 <- bplapply(res3fitstats, function(x) {
            if (length(x) > 0) {
                if (!is.na(x)) {
                    output <- mapply(function(anv, out, Singular, Conv) {
                        return(c(anv, out, Singular, Conv))
                    }, x$anv, x$fit, x$Singular, x$Conv, SIMPLIFY = FALSE)
                    return(output)
                } else {
                    return(NA)
                }
            } else {
                return(NA)
            }
        }, BPPARAM = SerialParam(progressbar = verbose))
        res3fitstats3 <- lapply(res3fitstats2, function(x) {
            if (length(x) > 0) {
                if (!is.na(x)) {
                    y <- do.call(rbind, x)
                }
            }
        })
        res3 <- do.call(rbind, res3fitstats3)
        res3 <- as.data.frame(res3)
    }

    if (verbose) {
        cat("Correcting P values", sep = "\n")
    }
    if (p.adjust.mode == "all") {
        if (method != "lme") {
            res3 <- bplapply(res3, function(x) {
                x$padj <- p.adjust(x$pval, method = p.adjust.method)
                row.names(x) <- x[, 1]
                x <- x[, -1, drop = FALSE]
                return(x)
            }, BPPARAM = SerialParam(progressbar = verbose))
        } else {
            p_cols <- grep("_P_", colnames(res3), value = TRUE)
            for (p in p_cols) {
                res3[, gsub("_P_", "_Padj_", p)] <- p.adjust(res3[, p], method = p.adjust.method)
            }
        }
    } else if (p.adjust.mode == "celltype") {
        if (method != "lme") {
            res3 <- bplapply(res3, function(x) {
                tmp <- split(x, x$celltype)
                tmp <- bplapply(tmp, function(y) {
                    y$padj <- p.adjust(y$pval, method = p.adjust.method)
                    row.names(y) <- y[, 1]
                    y <- y[, -1, drop = FALSE]
                    return(y)
                }, BPPARAM = SerialParam())
                tmp <- do.call(rbind, tmp)
                return(tmp)
            }, BPPARAM = SerialParam(progressbar = verbose))
        } else {
            celltypes <- split(res3, res3$celltypes)
            ctp <- bplapply(celltypes, function(x) {
                p_cols <- grep("_P_", colnames(x), value = TRUE)
                for (p in p_cols) {
                    x[, gsub("_P_", "_Padj_", p)] <- p.adjust(x[, p], method = p.adjust.method)
                }
                return(x)
            }, BPPARAM = SerialParam())
            res3 <- do.call(rbind, ctp)
        }
    }
    return(res3)
}
