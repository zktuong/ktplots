#' Plotting cellphonedb results
#'
#' @param cell_type1 cell type 1
#' @param cell_type2 cell type 2
#' @param scdata single-cell data. can be seurat/summarizedexperiment object
#' @param idents vector holding the idents for each cell or column name of scdata's metadata. MUST match cpdb's columns
#' @param means object holding means.txt from cpdb output
#' @param pvals object holding pvals.txt from cpdb output
#' @param max_size max size of points.
#' @param p.adjust.method correction method. p.adjust.methods of one of these options: c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none")
#' @param keep_significant_only logical. Default is FALSE. Switch to TRUE if you only want to plot the significant hits from cpdb.
#' @param split.by column name in the metadata/coldata table to split the spots by. Can only take columns with binary options. If specified, name to split by MUST be specified in the meta file provided to cpdb prior to analysis.
#' @param gene.family default = NULL. some predefined group of genes. can take one of these options: "chemokines", "Th1", "Th2", "Th17", "Treg", "costimulatory", "coinhibitory", "niche"
#' @param genes default = NULL. can specify custom list of genes if gene.family is NULL
#' @param scale logical. scale the expression to mean +/- SD. NULL defaults to TRUE.
#' @param standard_scale logical. scale the expression to range from 0 to 1. NULL defaults to FALSE.
#' @param col_option specify plotting colours
#' @param default_stlye default = TRUE. Show all mean values and trace significant interactions with `higlight` colour. If FALSE, significant interactions will be presented as a white ring.
#' @param noir default = FALSE. makes it b/w
#' @param highlight colour for highlighting p <0.05
#' @param highlight_size stroke size for highlight if p < 0.05. if NULL, scales to -log10(pval).
#' @param separator separator to use to split between celltypes. Unless otherwise specified, the separator will be `>@<`. Make sure the idents and split.by doesn't overlap with this.
#' @param special_character_search_pattern search pattern if the cell type names contains special character. NULL defaults to "/|:|\\?|\\*|\\+|[\\]|\\(|\\)".
#' @param ... passes arguments to grep for cell_type1 and cell_type2.
#' @return ggplot dot plot object of cellphone db output
#' @examples
#' \donttest{
#' data(kidneyimmune)
#' data(cpdb_output)
#' plot_cpdb("B cell", "CD4T cell", kidneyimmune, 'celltype', means, pvals, split.by = "Experiment", genes = c("CXCL13", "CD274", "CXCR5"))
#' plot_cpdb("B cell", "CD4T cell", kidneyimmune, 'celltype', means, pvals, split.by = "Experiment", gene.family = 'chemokines')
#' }
#' @import viridis
#' @import ggplot2
#' @import reshape2
#' @export

plot_cpdb <- function(cell_type1, cell_type2, scdata, idents, means, pvals, max_size = 8, p.adjust.method = NULL, keep_significant_only = FALSE, split.by = NULL, gene.family = NULL, genes = NULL, scale = NULL, standard_scale = NULL, col_option = viridis::viridis(50), default_style = TRUE, noir = FALSE, highlight = "red", highlight_size = NULL, separator = NULL, special_character_search_pattern = NULL, ...) {
    if (class(scdata) %in% c("SingleCellExperiment", "SummarizedExperiment")) {
        cat("data provided is a SingleCellExperiment/SummarizedExperiment object", sep = "\n")
        cat("extracting expression matrix", sep = "\n")
        requireNamespace('SummarizedExperiment')
        requireNamespace('SingleCellExperiment')
        exp_mat <- SummarizedExperiment::assay(scdata)
        metadata <- SummarizedExperiment::colData(scdata)
    } else if (class(scdata) == "Seurat") {
        requireNamespace('Seurat')
        cat("data provided is a Seurat object", sep = "\n")
        cat("extracting expression matrix", sep = "\n")
        exp_mat <- tryCatch(scdata@data, error = function(e) {
            tryCatch(Seurat::GetAssayData(object = scdata), error = function(e) {
                stop(paste0("are you sure that your data is normalized?"))
                return(NULL)
            })
        })
        metadata <- scdata@meta.data
    }

    if (length(separator) > 0){
        sep = separator
    } else {
        sep = '>@<'
    }

    means_mat <- means
    pvals_mat <- pvals
    rownames(means_mat) <- make.names(means_mat$interacting_pair, unique = TRUE)
    rownames(pvals_mat) <- make.names(pvals_mat$interacting_pair, unique = TRUE)
    colnames(means_mat) <- gsub("\\|", sep, colnames(means_mat))
    rownames(means_mat) <- gsub("_", "-", rownames(means_mat))
    rownames(means_mat) <- gsub("[.]", " ", rownames(means_mat))
    colnames(pvals_mat) <- gsub("\\|", sep, colnames(pvals_mat))
    rownames(pvals_mat) <- gsub("_", "-", rownames(pvals_mat))
    rownames(pvals_mat) <- gsub("[.]", " ", rownames(pvals_mat))

    if(length(p.adjust.method) > 0){
        pvals_tmp <- pvals[,12:ncol(pvals)]
        pvals_adj <- matrix(p.adjust(as.vector(as.matrix(pvals_tmp)), method=p.adjust.method),ncol=ncol(pvals_tmp))
        colnames(pvals_adj) <- colnames(pvals_tmp)
        pvals <- cbind(pvals[,c(1:11)], pvals_adj)
    }

    if (is.null(special_character_search_pattern)){
        pattern <- "/|:|\\?|\\*|\\+|[\\]|\\(|\\)|\\/"
    } else {
        pattern <- special_character_search_pattern
    }

    cell_type1_tmp <- unlist(strsplit(cell_type1, '*'))
    cell_type2_tmp <- unlist(strsplit(cell_type2, '*'))

    if (any(grepl(pattern, cell_type1_tmp))){
        idxz <- grep(pattern, cell_type1_tmp)
        cell_type1_tmp[idxz] <- paste0("\\", cell_type1_tmp[idxz])
        cell_type1 <- do.call(paste, c(as.list(cell_type1_tmp), sep = ""))
    } else {
        cell_type1 <- cell_type1
    }
    if (any(grepl(pattern, cell_type2_tmp))){
        idxz <- grep(pattern, cell_type2_tmp)
        cell_type2_tmp[idxz] <- paste0("\\", cell_type2_tmp[idxz])
        cell_type2 <- do.call(paste, c(as.list(cell_type2_tmp), sep = ""))
    } else {
        cell_type2 <- cell_type2
    }

    if(length(idents) > 1){
        ct1 = grep(cell_type1, idents, value = TRUE, ...)
        ct2 = grep(cell_type2, idents, value = TRUE, ...)
        checklabels1 <- any(idents %in% c(ct1,ct2))
    } else {
        ct1 = grep(cell_type1, metadata[[idents]], value = TRUE, ...)
        ct2 = grep(cell_type2, metadata[[idents]], value = TRUE, ...)
        checklabels1 <- any(metadata[[idents]] %in% c(ct1,ct2))
    }

    if(!is.null(split.by)){
        if(length(idents) > 1){
            labels <- paste0(metadata[[split.by]], "_", idents)
        } else {
            labels <- paste0(metadata[[split.by]], "_", metadata[[idents]])
        }

        labels <- factor(labels)
        labels <- levels(labels)
        groups <- factor(metadata[[split.by]])
        groups <- levels(groups)

        if(length(groups) > 0){
            # the purpose for this step is to allow for special characters to be used in the celltype grepping
            if(length(groups) > 1){
                labels2 = gsub(paste0(paste0(groups,'_'), collapse = '|'), '', labels)
            } else {
                labels2 = gsub(paste0(groups,'_'), '', labels)
            }

            # this returns the indices from the labels
            ct1 = grep(cell_type1, labels2, value = TRUE, ...)
            ct2 = grep(cell_type2, labels2, value = TRUE, ...)
        } else {
            if(length(idents) > 1){
                labels <- idents
            } else {
                labels <- metadata[[idents]]
            }
            labels <- factor(labels)
            labels <- levels(labels)

            ct1 = grep(cell_type1, labels, value = TRUE, ...)
            ct2 = grep(cell_type2, labels, value = TRUE, ...)
        }
    } else {
        if(length(idents) > 1){
            labels <- idents
        } else {
            labels <- metadata[[idents]]
        }
        labels <- factor(labels)
        labels <- levels(labels)

        ct1 = grep(cell_type1, labels, value = TRUE, ...)
        ct2 = grep(cell_type2, labels, value = TRUE, ...)
        ct1 = paste0(ct1, collapse = '|')
        ct2 = paste0(ct2, collapse = '|')
    }

    x1 = ct1[ct1 %in% '']
    x2 = ct2[ct2 %in% '']
    if (length(x1) > 0){
        ct1[ct1 %in% ''] <- NA
    }
    if (length(x2) > 0){
        ct2[ct2 %in% ''] <- NA
    }
    checklabels2 <- any(colnames(means_mat) %in% c(ct1,ct2))

    if(!checklabels1){
        if(length(idents) > 1){
            # relatively relaxed criteria to allow for the program to continue
            options(warn=-1)
            ct_1 <- grep(cell_type1, idents, value = TRUE, ...)
            ct_2 <- grep(cell_type2, idents, value = TRUE, ...)
            options(warn=0)
            checklabels2 <- any(idents %in% ct_1)
            if(checklabels2){
                checklabels2 <- any(idents %in% ct_2)
                if(!checklabels2){
                    stop('Cannot find cell types.\nThe error is mismatch between cell_type2 and the single cell metadata (or idents provided).')
                }
            } else {
                stop('Cannot find cell types.\nThe error is mismatch between cell_type1 and the single cell metadata (or idents provided).')
            }
        } else {
            options(warn=-1)
            ct_1 <- grep(cell_type1, metadata[[idents]], value = TRUE, ...)
            ct_2 <- grep(cell_type2, metadata[[idents]], value = TRUE, ...)
            options(warn=0)
            checklabels2 <- any(metadata[[idents]] %in% ct_1)
            if(checklabels2){
                checklabels2 <- any(metadata[[idents]] %in% ct_2)
                if(!checklabels2){
                    stop('Cannot find cell types.\nThe error is mismatch between cell_type2 and the single cell metadata (or idents provided).')
                }
            } else {
                stop('Cannot find cell types.\nThe error is mismatch between cell_type1 and the single cell metadata (or idents provided).')
            }
        }
    }

    if(!checklabels2){
        # relatively relaxed criteria to allow for the program to continue
        options(warn=-1)
        ct_1 <- grep(cell_type1, colnames(means_mat), value = TRUE)
        ct_2 <- grep(cell_type2, colnames(means_mat), value = TRUE)
        options(warn=0)
        checklabels2 <- any(colnames(means_mat) %in% ct_1)
        if(checklabels2){
            checklabels2 <- any(colnames(means_mat) %in% ct_2)
            if(!checklabels2){
                stop('Cannot find cell types. The error is mismatch between cell_type2 and the cpdb metadata.')
            }
        } else {
            stop('Cannot find cell types. The error is mismatch between cell_type1 and the cpdb metadata.')
        }
    }

    if(checklabels1 & checklabels2){
        cat('Found cell types in the input data provided. Proceeding with plotting.', sep = "\n")
    }

    if (is.null(gene.family) & is.null(genes)){
        cat("options genes or gene.family are not specified.\nusing entire cpdb output.\n")
        query <- grep('', means_mat$interacting_pair)
        cat("for future reference, genes or gene.family can be specified, not both.\ngene.family can be one of the following:", sep = "\n")
        print(c("chemokines", "Th1", "Th2", "Th17", "Treg", "costimulatory", "coinhibitory", "niche"))
        cat("otherwise, please provide gene(s) as a vector in the genes option", sep = "\n")
    }

    if (!is.null(gene.family) & !is.null(genes)){
        stop("Please specify either genes or gene.family, not both")
        cat("gene.family can be one of the following:", sep = "\n")
        print(c("chemokines", "Th1", "Th2", "Th17", "Treg", "costimulatory", "coinhibitory", "niche"))
        cat("otherwise, please provide gene(s) as a vector in the genes option", sep = "\n")
    }

    if(!is.null(gene.family) & is.null(genes)){
        chemokines <- grep("^CXC|CCL|CCR|CX3|XCL|XCR", means_mat$interacting_pair)
        th1 <- grep("IL2|IL12|IL18|IL27|IFNG|IL10|TNF$|TNF |LTA|LTB|STAT1|CCR5|CXCR3|IL12RB1|IFNGR1|TBX21|STAT4", means_mat$interacting_pair)
        th2 <- grep("IL4|IL5|IL25|IL10|IL13|AREG|STAT6|GATA3|IL4R", means_mat$interacting_pair)
        th17 <- grep("IL21|IL22|IL24|IL26|IL17A|IL17A|IL17F|IL17RA|IL10|RORC|RORA|STAT3|CCR4|CCR6|IL23RA|TGFB", means_mat$interacting_pair)
        treg <- grep("IL35|IL10|FOXP3|IL2RA|TGFB", means_mat$interacting_pair)
        costimulatory <- grep("CD86|CD80|CD48|LILRB2|LILRB4|TNF|CD2|ICAM|SLAM|LT[AB]|NECTIN2|CD40|CD70|CD27|CD28|CD58|TSLP|PVR|CD44|CD55|CD[1-9]", means_mat$interacting_pair)
        coinhibitory <- grep("SIRP|CD47|ICOS|TIGIT|CTLA4|PDCD1|CD274|LAG3|HAVCR|VSIR", means_mat$interacting_pair)
        niche <- grep("CSF", means_mat$interacting_pair)
        query_group <- list(chemokines = chemokines, chemokine = chemokines, th1 = th1, th2 = th2, th17 = th17, treg = treg, costimulatory = costimulatory, coinhibitory = coinhibitory, costimulation = costimulatory, coinhibition = coinhibitory, niche = niche)
    } else if (is.null(gene.family) & !is.null(genes)){
        query <- grep(paste(genes, collapse="|"), means_mat$interacting_pair)
    }

    create_celltype_query <- function(ctype1, ctype2, sep){
        ct1 = list()
        ct2 = list()
        for (i in 1:length(ctype2)){
            ct1[i] = paste0('^', ctype1, sep, ctype2[i], '$')
            ct2[i] = paste0('^', ctype2[i], sep, ctype1, '$')
        }
        ct_1 = do.call(paste0, list(ct1, collapse = '|'))
        ct_2 = do.call(paste0, list(ct2, collapse = '|'))
        ct = list(ct_1, ct_2)
        ct = do.call(paste0, list(ct, collapse = '|'))
        return(ct)
    }

    keep_interested_groups <- function(g, ct, sep){
        ctx <- strsplit(ct, "\\|")[[1]]
        idx <- grep(paste0(g, paste0(".*",sep), g), ctx)
        ctx <- ctx[idx]
        ctx <- paste0(ctx, collapse = "|")
        return(ctx)
    }

    if(!is.null(split.by)){
        if(length(idents) > 1){
            labels <- paste0(metadata[[split.by]], "_", idents)
        } else {
            labels <- paste0(metadata[[split.by]], "_", metadata[[idents]])
        }

        chk1 = class(metadata[[split.by]])
        chk2 = class(metadata[[idents]])
        if (chk1 == 'factor' & chk2 == 'factor'){
            labels <- factor(labels, levels = paste0(levels(metadata[[split.by]]), '_', rep(levels(metadata[[idents]]), each = length(levels(metadata[[split.by]])))))
        } else {
            labels <- factor(labels)
        }
        labels <- levels(labels)
        groups <- factor(metadata[[split.by]])
        groups <- levels(groups)

        if(length(groups) > 0){
                # the purpose for this step is to allow for special characters to be used in the celltype grepping
            if(length(groups) > 1){
                labels2 = gsub(paste0(paste0(groups,'_'), collapse = '|'), '', labels)
            } else {
                labels2 = gsub(paste0(groups,'_'), '', labels)
            }
                # this returns the indices from the labels
            ct1 = grep(cell_type1, labels2, ...)
            ct2 = grep(cell_type2, labels2, ...)

            c_type1 <- as.list(labels[ct1])
            c_type2 <- as.list(labels[ct2])

            grp <- as.list(groups)
            celltype = list()
            for (i in 1:length(c_type1)){
                celltype[[i]] <- create_celltype_query(c_type1[[i]], c_type2, sep)
                celltype[[i]] <- lapply(grp, keep_interested_groups, celltype[[i]], sep)
            }

            for (i in 1:length(celltype)){
                celltype[[i]] <- celltype[[i]][-which(celltype[[i]] == '')]
            }

            celltype <- lapply(celltype, unlist)
            if (any(unlist(lapply(celltype, is.null)))){
                rm <- which(unlist(lapply(celltype, is.null)))
                celltype <- celltype[-rm]
            }
            celltype2 <- lapply(celltype, function(x) {
                xx <- as.list(unlist(strsplit(x, "\\|")))
                xx <- lapply(xx, function(z) {
                    xz <- paste0('^', z, '$')
                    tmp_ <- unlist(strsplit(xz, '*'))
                    if (any(grepl(pattern, tmp_))){
                        idxz <- grep(pattern, tmp_)
                        tmp_[idxz] <- paste0("\\", tmp_[idxz])
                    }
                    xz <- do.call(paste, c(as.list(tmp_), sep = ""))
                    return(xz)})
                zz <- paste0(unlist(xx), collapse = "|")
                return(zz)
            })

            cell_type <- do.call(paste0, list(celltype2, collapse = "|"))

        } else {
            if(length(idents) > 1){
                labels <- idents
            } else {
                labels <- metadata[[idents]]
            }
            labels <- factor(labels)
            labels <- levels(labels)
            c_type1 = as.list(grep(cell_type1, labels, value = TRUE, ...))
            c_type2 = as.list(grep(cell_type2, labels, value = TRUE, ...))

            celltype = list()
            for (i in 1:length(c_type1)){
                celltype[[i]] <- create_celltype_query(c_type1[[i]], c_type2, sep)
            }
            cell_type <- do.call(paste0, list(celltype, collapse = "|"))
        }
    } else {
        if(length(idents) > 1){
            labels <- idents
        } else {
            labels <- metadata[[idents]]
        }
        labels <- factor(labels)
        labels <- levels(labels)
        c_type1 = as.list(grep(cell_type1, labels, value = TRUE))
        c_type2 = as.list(grep(cell_type2, labels, value = TRUE))

        celltype = list()
        for (i in 1:length(c_type1)){
            celltype[[i]] <- create_celltype_query(c_type1[[i]], c_type2, sep)
        }
        cell_type <- do.call(paste0, list(celltype, collapse = "|"))
        cell_type_tmp <- unlist(strsplit(cell_type, '*'))
        if (any(grepl(pattern, cell_type_tmp))){
            idxz <- grep(pattern, cell_type_tmp)
            cell_type_tmp[idxz] <- paste0("\\", cell_type_tmp[idxz])
            cell_type <- do.call(paste, c(as.list(cell_type_tmp), sep = ""))
        } else {
            cell_type <- cell_type
        }
    }


    if(!is.null(gene.family) & is.null(genes)){
        means_mat <- means_mat[query_group[[tolower(gene.family)]], grep(cell_type, colnames(means_mat), ...)]
        pvals_mat <- pvals_mat[query_group[[tolower(gene.family)]], grep(cell_type, colnames(pvals_mat), ...)]
    } else if (is.null(gene.family) & !is.null(genes) | is.null(gene.family) & is.null(genes)){
        means_mat <- means_mat[query, grep(cell_type, colnames(means_mat), ...)]
        pvals_mat <- pvals_mat[query, grep(cell_type, colnames(pvals_mat), ...)]
    }

    if (length(means_mat) == 0){
        stop('Please check your options for split.by and your celltypes.')
    }

    # rearrange the columns so that it interleaves the two groups
    if(!is.null(split.by)){
        if(length(groups) > 0){
            grp <- as.list(groups)
            group_i <- lapply(grp, function(g){
                gx <- grep(g, colnames(means_mat), ...)
                return(gx)
            })
            group_id <- do.call(rbind, group_i)
            means_mat <- means_mat[,as.vector(group_id)]
            if (dim(pvals_mat)[2] > 0){
                pvals_mat <- pvals_mat[,as.vector(group_id)]
            } else {
                stop('No significant hits.')
            }
        }
    }

    if (keep_significant_only){
        if (dim(pvals_mat)[2] == 0){
            stop('No significant hits.')
        }
    }

    if(nrow(means_mat) > 2){
        d <- dist(as.data.frame(means_mat))
        h <- hclust(d)
        means_mat <- means_mat[h$order, ]
        pvals_mat <- pvals_mat[h$order, ]
    }

    # scaling
    if (length(standard_scale) > 0){
        if (standard_scale){
            means_mat_ <- apply(means_mat,1,range01)
            means_mat_ <- t(means_mat_)
        } else {
            means_mat_ <- means_mat
        }
    }

    if(length(scale) < 1){
        if(length(standard_scale) > 0){
            if (standard_scale){
                means_mat2 <- means_mat_
            } else {
                means_mat2 <- means_mat
            }
        } else {
            means_mat2 <- t(scale(t(means_mat)))
        }
    } else {
        if (scale){
            if(length(standard_scale) > 0){
                if (standard_scale){
                    means_mat2 <- means_mat_
                } else {
                    means_mat2 <- t(scale(t(means_mat)))
                }
            } else {
                means_mat2 <- t(scale(t(means_mat)))
            }
        } else {
            means_mat2 <- means_mat
        }
    }

    pvals_mat2 <- as.matrix(pvals_mat)
    means_mat2 <- as.matrix(means_mat2)
    xx <- which(means_mat == 0)
    if (length(xx) > 0){
        means_mat2[which(means_mat == 0)] <- NA
    }

    # remove rows that are entirely NA
    pvals_mat2 <- pvals_mat2[rowSums(is.na(means_mat2)) != ncol(means_mat2), ,drop = FALSE]
    means_mat2 <- means_mat2[rowSums(is.na(means_mat2)) != ncol(means_mat2), ,drop = FALSE]

    if ((length(standard_scale) > 0 && standard_scale) | (length(scale) > 0 && scale) | (length(scale) < 1 && length(standard_scale) < 1)){
        df_means <- melt(means_mat2, value.name = "scaled_means")
    } else {
        df_means <- melt(means_mat2, value.name = "means")
    }

    if(length(p.adjust.method) > 0){
        df_pvals <- melt(pvals_mat2, value.name = "pvals_adj")
        df <- data.frame(cbind(df_means, pvals_adj = df_pvals$pvals_adj))
        xp <- which(df$pvals_adj == 1)
        if (length(xp) > 0){
            df$pvals_adj[which(df$pvals_adj == 1)] <- NA
        }
        if (keep_significant_only){
            # keep the entire row/ all the comparisons
            df_ <- split(df, as.character(df$Var1))
            anysig <- lapply(df_, function(x){
                keep <- any(x$pvals_adj < 0.05)
                return(keep)
            })
            df_ <- df_[which(unlist(anysig))]
            names(df_) <- NULL
            df <- do.call(rbind, df_)
        }
        df$pvals_adj[which(df$pvals_adj == 0)] <- 0.001
    } else {
        df_pvals <- melt(pvals_mat2, value.name = "pvals")
        df <- data.frame(cbind(df_means, pvals = df_pvals$pvals))
        xp <- which(df$pvals == 1)
        if (length(xp) > 0){
            df$pvals[which(df$pvals == 1)] <- NA
        }
        if (keep_significant_only){
            # keep the entire row/ all the comparisons
            df_ <- split(df, as.character(df$Var1))
            anysig <- lapply(df_, function(x){
                keep <- any(x$pvals < 0.05)
                return(keep)
            })
            df_ <- df_[which(unlist(anysig))]
            names(df_) <- NULL
            df <- do.call(rbind, df_)
        }
        df$pvals[which(df$pvals == 0)] <- 0.001
    }

    if(!is.null(split.by)){
        if(length(groups) > 0){
            grp <- as.list(groups)
            grp2 <- lapply(grp, function(i){
                x <- paste0(i,'_')
                return(x)
            })
            searchterm <- do.call(paste, list(grp2, collapse = "|"))
            df$group <- gsub(searchterm, "", df$Var2)
        } 
    } else {
        df$group <- df$Var2
    }

    if (keep_significant_only){
        if ((length(standard_scale) > 0 && standard_scale) | (length(scale) > 0 && scale) | (length(scale) < 1 && length(standard_scale) < 1)){
            if (length(df$scaled_means) == 0){
                stop('No significant genes found and plotting will not proceed.')
            }
        } else {
            if (length(df$means) == 0){
                stop('No significant genes found and plotting will not proceed.')
            }
        }
    }

    df$Var2 <- gsub(sep, '-', df$Var2)
    final_levels = unique(df$Var2)
    df$Var2 <- factor(df$Var2, unique(df$Var2))
    df$x_means_ <- df[,colnames(df_means)[3]]
    df$x_means_[df[,colnames(df)[4]] < 0.05] <- NA
    df$x_stroke = df$x_means_
    df$x_stroke[!is.na(df$x_stroke)] <- 0
    df$x_stroke[is.na(df$x_stroke)] <- 2

    if (default_style){
        if ((length(standard_scale) > 0 && standard_scale) | (length(scale) > 0 && scale) | (length(scale) < 1 && length(standard_scale) < 1)){
            if(length(p.adjust.method) > 0 && p.adjust.method != 'none'){
                g <- ggplot(df, aes(x = Var2, y = Var1, color = -log10(pvals_adj), fill = scaled_means, size = scaled_means))
            } else {
                g <- ggplot(df, aes(x = Var2, y = Var1, color = -log10(pvals), fill = scaled_means, size = scaled_means))
            }
        } else {
            if(length(p.adjust.method) > 0 && p.adjust.method != 'none'){
                g <- ggplot(df, aes(x = Var2, y = Var1, color = -log10(pvals_adj), fill = means, size = means))
            } else {
                g <- ggplot(df, aes(x = Var2, y = Var1, color = -log10(pvals), fill = means, size = means))
            }
        }

        if (!is.null(highlight_size)){
            g <- g + geom_point(pch = 21, na.rm=TRUE, stroke = highlight_size)
        } else {
            if(length(p.adjust.method) > 0 && p.adjust.method != 'none'){
                s = -log10(df$pvals_adj)
                s[is.na(s)] <- 0
                g <- g + geom_point(pch = 21, na.rm=TRUE, stroke = s)
            } else {
                s = -log10(df$pvals)
                s[is.na(s)] <- 0
                g <- g + geom_point(pch = 21, na.rm=TRUE, stroke = s)
            }
        }
        g <- g +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 45, hjust = 0, color = '#000000'),
            axis.text.y = element_text(color = '#000000'),
            axis.ticks = element_blank(),
            axis.title.x = element_blank(),
            axis.title.y = element_blank()) +
        scale_x_discrete(position = "top") +
        scale_color_gradientn(colors = highlight, na.value = "white") +
        scale_radius(range = c(0,max_size))

        if(noir){
            g <- g + scale_fill_gradient(low = "white", high = "#131313", na.value = "white")
        } else {
            if(length(col_option) == 1){
                g <- g + scale_fill_gradientn(colors = colorRampPalette(c("white", col_option))(100), na.value = "white")
            } else {
                g <- g + scale_fill_gradientn(colors = c("white", colorRampPalette(col_option)(99)), na.value = "white")
            }
        }
    } else {
        if ((length(standard_scale) > 0 && standard_scale) | (length(scale) > 0 && scale) | (length(scale) < 1 && length(standard_scale) < 1)){
            g <- ggplot(df, aes(x = Var2, y = Var1, size = scaled_means, color = scaled_means))
            df2 <- df %>% filter(is.na(x_means_))
            g <- g + geom_point(pch = 16, na.rm = TRUE)
            g <- g + geom_point(data = df2, aes(x = Var2, y = Var1, size = scaled_means, fill = x_means_, stroke = x_stroke), pch = 21, na.rm = TRUE)
        } else {
            g <- ggplot(df, aes(x = Var2, y = Var1, size = means, color = means))
            df2 <- df %>% filter(is.na(x_means_))
            g <- g + geom_point(pch = 16, na.rm = TRUE)
            g <- g + geom_point(data = df2, aes(x = Var2, y = Var1, size = scaled_means, fill = x_means_, stroke = x_stroke), pch = 21, na.rm = TRUE)
        }
        g <- g +
        theme_bw() +
        scale_fill_gradientn(colors = col_option, na.value = 'white', guide = FALSE) + scale_colour_gradientn(colors = col_option) +
        theme(axis.text.x = element_text(angle = 45, hjust = 0, color = '#000000'),
            axis.text.y = element_text(color = '#000000'),
            axis.ticks = element_blank(),
            axis.title.x = element_blank(),
            axis.title.y = element_blank()) +
        scale_x_discrete(position = "top") +
        scale_radius(range = c(0,max_size))
    }

    if(!is.null(gene.family) & is.null(genes)){
        g <- g + ggtitle(gene.family)
    }
    return(g)
}
