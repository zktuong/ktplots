#' Plotting cellphonedb results
#'
#' @param meta data.frame containing the sample name, cellphonedb 'out' folder path (containing means.txt and pvalues.txt), and single-cell object file path (.h5ad or .rds)
#' @param celltypes subset celltypes for comparison
#' @param method one of 'ttest', 'wilcox', 'lme'
#' @param p.adjust.method defaults to p.adjust methods
#' @param groupby for significance testing. only if method is t.test or wilcoxon.
#' @param formula for signfiicance testing. only if method is lme.
#' @param BPPARAM BiocParallelParam class.
#' @param verbose Whether or not to print messages.
#' @param ... passes arguments to significance tests (pairwise.t.test, pairwise.wilcox.text, lmerTest::lmer)
#' @return results for plotting
#' @examples
#' \donttest{}
#' @import BiocParallel
#' @import dplyr
#' @export
compare_cpdb <- function(meta, celltypes, method = c('ttest', 'wilcox', 'lme'), p.adjust.method = 'fdr', groupby = NULL, formula = NULL, BPPARAM = SerialParam(), verbose = TRUE) {
    options(warn = -1)
    sample <- meta[, 1]
    cpdb_out_folder <- meta[, 2]
    expression_file <- meta[, 3]
    names(cpdb_out_folder) <- sample
    names(expression_file) <- sample

    if (verbose) {
        cat("Reading cellphonedb outputs", sep = "\n")
    }
    means <- bplapply(cpdb_out_folder, function(x) read.delim(paste0(x, "means.txt"),
        check.names = FALSE), BPPARAM = SerialParam(progressbar = verbose))
    pvals <- bplapply(cpdb_out_folder, function(x) read.delim(paste0(path, x, "pvalues.txt"),
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
    if (verbose){
        cat(paste0("Running for ", dim(combs)[1], " celltype combinations"), sep = "\n")
    }

    get_means <- function(x){
    if (!is.na(x)){
        out <- x$means
        names(out) <- paste0(x$group, '>@<', x$Var1)
        if (any(is.na(out))){
            out <- out[!is.na(out)]
        }
        return(out)
    } else {
        return(NA)
    }
    }

    # extract interactions with plot_cpdb
    if (verbose){
        cat('Extracting cpdb results across combinations', sep = "\n")
        start_time <- Sys.time()
    }

    ct_sigs <- bplapply(comb_list, function(x){
        bpmapply(function(sce, mean, pval) {
                s <- tryCatch(plot_cpdb(cell_type1 = paste0(x[1], '$'),
                    cell_type2 = paste0(x[2], '$'),
                    scdata = sce,
                    idents = idents, # column name where the cell ids are located in the metadata
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
    if (verbose){
        end_time <- Sys.time()
        end_time - start_time
    }

    # collapse to matrix
    requireNamespace("plyr")
    meta <- bplapply(ct_sigs, function(y) {
        z <- as.data.frame(plyr::rbind.fill(lapply(y, function(y) {
            tryCatch(as.data.frame(t(y), stringsAsFactors = FALSE), error = function(e) return(NA))
        })))
        row.names(z) <- names(sces)
        return(z)
    }, BPPARAM = BPPARAM)

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
    }, BPPARAM = SerialParam())

    # split each column of into a list
    res2 <- bplapply(meta2, function(x) {
        tmp <- as.list(x)
        tmp <- bplapply(tmp, function(z) {
            z <- as.numeric(z)
            names(z) <- row.names(x)
            return(z)
        }, BPPARAM = BPPARAM)
        return(tmp)
    }, BPPARAM = SerialParam())

    requireNamespace('reshape2')
    test_fun <- function(int_score, data, col, method = c("t.test", "wilcox"), ...) {
        # 'C'
        data$int_score <- int_score
        if (class(data[, col]) == "factor") {
            lvls <- levels(data[, col])
        } else {
            data[, col] <- factor(data[, col])
            lvls <- levels(data[, col])
        }

        if (method == "wilcox") {
            test <- pairwise.wilcox.test(data$int_score, data[, col], p.adjust.method = "none",
                ...)
        } else if (method == "t.test") {
            test <- pairwise.t.test(data$int_score, data[, col], p.adjust.method = "none",
                ...)
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

    if (verbose){
        cat('Running significance testing', sep = "\n")
        start_time <- Sys.time()
    }

    if (method %in% c('t.test', 'wilcox')){
        if (is.null(groupby)) {
            stop("Please provide column name for contrasts to be extracted.")
        }
        res3 <- bplapply(res2, function(x) {
            tmp <- bplapply(x, function(y) {
                suppressMessages(suppressWarnings(tryCatch(test_fun(int_score = y,
                    data = metadata, col = groupby, method = method, ...), error = function(e) return(NA))))
            }, BPPARAM = BPPARAM)
            tmp <- plyr::ldply(tmp, data.frame)
            return(tmp)
        }, BPPARAM = SerialParam(progressbar = verbose))
        res3 <- do.call(rbind, res3)
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
            cat("running lmer model", sep = "\n")
        }
        res3_ <- bplapply(res2, function(x) {
            bplapply(x, function(int_score) {
                suppressMessages(suppressWarnings(tryCatch(lmer(fullFormula, data = metadata,
                    ...), error = function(e) return(NA))))
            }, BPPARAM = BPPARAM)
        }, BPPARAM = SerialParam(progressbar = verbose))

        if (verbose) {
            cat("extracting statistics", sep = "\n")
        }

        res3fitstats <- bplapply(res3_, function(fit) {
            if (length(fit) > 0) {
                if (!is.na(fit)) {
                    out <- lapply(fit, function(fit2) {
                        if (!is.na(fit2)) {
                          y <- suppressWarnings(summary(fit2))
                          E <- y$coefficients[-1, 1]
                          names(E) <- paste0("fit_estimates_", names(E))
                          P <- y$coefficients[-1, 5]
                          names(P) <- paste0("fit_P_", names(P))
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
                    wald <- lapply(fit, function(fit2) {
                        if (!is.na(fit2)) {
                          anova <- car::Anova(fit2)
                          wp <- anova[, 3]
                          names(wp) <- paste0("Wald_P_", row.names(anova))
                          return(wp)
                        } else {
                          return(NA)
                        }
                    })
                    return(list(wald = wald, fit = out, Singular = Singular, Conv = Conv))
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
                    output <- mapply(function(wald, out, Singular, Conv) {
                        return(c(wald, out, Singular, Conv))
                    }, x$wald, x$fit, x$Singular, x$Conv, SIMPLIFY = FALSE)
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
        res3 <- do.call(rbind, res3s2)

    }
    if (verbose){
        end_time <- Sys.time()
        end_time - start_time
    }

    if (p.adjust.method != 'none'){
        if (verbose){
            cat('correcting P values', sep = '\n')
        }
        if (method %in% c('t.test', 'wilcox')){
            res3 <- bplapply(res3, function(x) {
                        x$padj <- p.adjust(x$pval, method = p.adjust.method)
                        row.names(x) <- x[,1]
                        x <- x[,-1]
                        return(x)}, BPPARAM = SerialParam(progressbar = verbose))
        } else {
            p_cols <- grep('_P_', colnames(res3), value = TRUE)
            for (p in p_cols)){
                res3[,gsub('_P_', '_Q_', i)] <- p.adjust(res3[,p], method = p.adjust.method)
            }
        }

    }

    return(res3)
}
