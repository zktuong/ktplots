#' Plotting correlations on Vissium data in Seurat
#'
## main function
#' @name correlationSpot
#' @param st spatial transcriptomics data in Seurat.
#' @param genes gene or genes of interest for performing correlations. Must exist as row name(s) in the `rna_slot`.
#' @param celltypes celltype or celltypes of interest for performing correlations. Must exist as row name(s) in the `label_slot`.
#' @param geneset geneset or column in meta.data for performing correlations. Must exist as column name(s) in the meta.data.
#' @param mode whether or not to restrict the output to just high expressing values, low expressing values or both. Caveat with using both is low expressing ~ low expressing will still return a high correlation value.
#' @param cutoff percentile cut off for determining mode output.
#' @param standardize whether or not to scale values from 0 to 1 before performing correlations.
#' @param dims number of dimensions used for calculating neighborhoods in spatial data.
#' @param k.params number of k neighbors for calculating neighborhoods in spatial data.
#' @param resolution resolution of spatial clustering.
#' @param rna_slot name of gene expression slot. Defaults to 'SCT'.
#' @param label_slot name of label/prediction slot. Defaults to 'predictions'.
#' @param by whether or not to define spatial clusters based on image spatial location or gene expression.
#' @param average_by_cluster whether or not to return the output averaged across clusters.
#' @param ... passed to Seurat::SpatialFeaturePlot
#' @return SpatialFeaturePlot
#' @export

correlationSpot <- function(
    st, genes = NULL, celltypes = NULL, geneset = NULL, mode = c(
        "high",
        "low", "both"
    ), cutoff = 0.5, standardize = TRUE, dims = 1:30, k.params = 10,
    resolution = 1, rna_slot = "SCT", label_slot = "predictions", by = c(
        "image",
        "expression"
    ), average_by_cluster = FALSE, ...) {
    requireNamespace("Seurat")
    requireNamespace("FNN")

    st1 <- st
    Seurat::DefaultAssay(st1) <- rna_slot
    st1 <- Seurat::FindNeighbors(st1, reduction = "pca", dims = dims, k.param = k.params)
    if (average_by_cluster) {
        st1 <- Seurat::FindClusters(st1, verbose = FALSE, resolution = resolution)
        st1 <- Seurat::AddMetaData(st1, Seurat::Idents(st1), "main_cluster")
    }

    # extract SNN graph
    if (by == "expression") {
        graph <- st1@graphs[[paste0(rna_slot, "_nn")]]
        knn_idx <- FNN::get.knn(graph, k = k.params)$nn.index
    } else if (by == "image") {
        mat <- st1@images$slice1@coordinates[, c("row", "col")]
        requireNamespace("stats")
        graph <- as.matrix(stats::dist(mat))
        knn_idx <- FNN::get.knn(graph, k = k.params)$nn.index
    }

    # proceed and extract the values
    if (!is.null(genes)) {
        rna_exp <- Seurat::GetAssayData(st1)[genes, ]
    } else {
        rna_exp <- NULL
        if (!is.null(geneset)) {
            genes <- geneset
            rna_exp <- t(st1@meta.data[, genes])
        } else {
            rna_exp <- NULL
        }
    }
    if (!is.null(celltypes)) {
        cellspot_pred <- st[[label_slot]][celltypes, ]
    } else {
        cellspot_pred <- NULL
    }

    if (!is.null(rna_exp) && !is.null(cellspot_pred)) {
        if (length(genes) > 1 && length(celltypes) > 1) {
            df <- t(rbind(as.matrix(rna_exp), as.matrix(cellspot_pred)))
        } else {
            if (length(celltypes) > 1 && length(genes) == 1) {
                df <- cbind(rna_exp, t(as.matrix(cellspot_pred)))
            } else if (length(celltypes) == 1 && length(genes) > 1) {
                df <- cbind(t(as.matrix(rna_exp)), as.numeric(cellspot_pred))
            } else {
                df <- data.frame(as.numeric(rna_exp), as.numeric(cellspot_pred))
            }
        }
        colnames(df) <- c(genes, celltypes)
    } else if (!is.null(cellspot_pred) && is.null(rna_exp)) {
        if (length(celltypes) > 1) {
            df <- t(as.matrix(cellspot_pred))
        } else {
            df <- data.frame(as.numeric(cellspot_pred))
        }
        colnames(df) <- celltypes
    } else if (!is.null(rna_exp) && is.null(cellspot_pred)) {
        if (length(genes) > 1) {
            df <- t(as.matrix(rna_exp))
        } else {
            df <- data.frame(as.numeric(rna_exp))
        }
        colnames(df) <- genes
    }

    if (standardize) {
        df <- apply(df, 2, range01)
    }

    # get the nth percentile of each column
    requireNamespace("stats")
    r1 <- apply(df, 2, stats::quantile, cutoff)
    if (mode != "both") {
        if (mode == "high") {
            for (i in seq_along(r1)) {
                df[df[, i] < r1[i], i] <- 0
            }
        } else if (mode == "low") {
            for (i in seq_along(r1)) {
                df[df[, i] > r1[i], i] <- 0
            }
        }
    }

    if (average_by_cluster) {
        if (ncol(df) == 2) {
            df <- cbind(as.data.frame(df), main_cluster = st1@meta.data[, c("main_cluster")])
        } else {
            requireNamespace("stats")
            if (length(genes) > 1 && length(celltypes) > 1) {
                dx_g <- stats::prcomp(df[, genes])
                dx_c <- stats::prcomp(df[, celltypes])
                dx_g <- as.data.frame(dx_g$x)[, 1]
                dx_c <- as.data.frame(dx_c$x)[, 1]
                df <- cbind(dx_g, dx_c, main_cluster = st1@meta.data[, "main_cluster"])
            } else if (length(genes) > 1 && length(celltypes) < 1) {
                dx_g <- stats::prcomp(df[, genes])
                dx_g <- as.data.frame(dx_g$x)[, 1]
                df <- cbind(dx_g, main_cluster = st1@meta.data[, "main_cluster"])
            } else if (length(celltypes) > 1 && length(genes) == 1) {
                dx_g <- df[, genes]
                dx_c <- stats::prcomp(df[, celltypes])
                dx_c <- as.data.frame(dx_c$x)[, 1]
                df <- cbind(dx_g, dx_c, main_cluster = st1@meta.data[, c("main_cluster")])
            } else if (length(celltypes) == 1 && length(genes) > 1) {
                dx_c <- df[, celltypes]
                dx_g <- stats::prcomp(df[, genes])
                dx_g <- as.data.frame(dx_g$x)[, 1]
                df <- cbind(dx_g, dx_c, main_cluster = st1@meta.data[, c("main_cluster")])
            } else if (length(celltypes) > 1 && length(genes) < 1) {
                dx_c <- stats::prcomp(df[, celltypes])
                dx_c <- as.data.frame(dx_c$x)[, 1]
                df <- cbind(dx_c, main_cluster = st1@meta.data[, "main_cluster"])
            }
        }
    } else {
        if (ncol(df) == 2) {
            df <- as.data.frame(df)
        } else {
            requireNamespace("stats")
            if (length(genes) > 1 && length(celltypes) > 1) {
                dx_g <- stats::prcomp(df[, genes])
                dx_c <- stats::prcomp(df[, celltypes])
                dx_g <- as.data.frame(dx_g$x)[, 1]
                dx_c <- as.data.frame(dx_c$x)[, 1]
                df <- cbind(dx_g, dx_c)
            } else if (length(celltypes) > 1 && length(genes) == 1) {
                dx_g <- df[, genes]
                dx_c <- stats::prcomp(df[, celltypes])
                dx_c <- as.data.frame(dx_c$x)[, 1]
                df <- cbind(dx_g, dx_c)
            } else if (length(celltypes) == 1 && length(genes) > 1) {
                dx_c <- dx_g[, celltypes]
                dx_g <- stats::prcomp(df[, genes])
                dx_g <- as.data.frame(dx_g$x)[, 1]
                df <- cbind(dx_g, dx_c)
            } else if (length(genes) > 1 && length(celltypes) < 1) {
                dx_g <- stats::prcomp(df[, genes])
                dx_g <- as.data.frame(dx_g$x)[, 1]
                df <- dx_g
            } else if (length(celltypes) > 1 && length(genes) < 1) {
                dx_c <- stats::prcomp(df[, celltypes])
                dx_c <- as.data.frame(dx_c$x)[, 1]
                df <- dx_c
            }
        }
    }

    dfl <- as.data.frame(df)
    nn_list <- list()
    for (i in seq_along(1:nrow(knn_idx))) {
        i_vec <- knn_idx[i, ]
        nn_list[[i]] <- dfl[i_vec, ]
    }

    if (average_by_cluster) {
        nn_list <- lapply(nn_list, function(x) {
            requireNamespace("stats")
            tmp <- x[, -which(colnames(x) %in% c("main_cluster"))]
            res <- stats::cor(tmp[, 1], tmp[, 2])
            return(res)
        })
    } else {
        nn_list <- lapply(nn_list, function(tmp) {
            requireNamespace("stats")
            res <- stats::cor(tmp[, 1], tmp[, 2])
            return(res)
        })
    }

    rv <- unlist(nn_list)
    names(rv) <- row.names(st1@meta.data)
    if (average_by_cluster) {
        dfl$cor <- rv
        dfl <- split(dfl, dfl$main_cluster)
        dfl <- lapply(dfl, function(x) mean(x$cor))
        rdf <- unlist(dfl)
        st1$correlation <- as.numeric(as.character(factor(st1$main_cluster, labels = rdf)))
        st$correlation <- st1$correlation
    } else {
        st <- Seurat::AddMetaData(st, rv, "correlation")
    }

    g <- Seurat::SpatialFeaturePlot(st, features = "correlation", ...)
    return(g)
}
