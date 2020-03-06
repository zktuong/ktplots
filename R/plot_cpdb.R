#' Plotting cellphonedb results
#'
#' @param cell_type1 cell type 1
#' @param cell_type2 cell type 2
#' @param scdata single-cell data. can be seurat/summarizedexperiment object
#' @param idents vector holding the idents for each cell or column name of scdata's metadata. MUST match cpdb's columns
#' @param means_file object holding means.txt from cpdb output
#' @param pvals_file object holding pvals.txt from cpdb output
#' @param p.adjust.method correction method. p.adjust.methods of one of these options: c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none")
#' @param keep_significant_only logical. Default is FALSE. Switch to TRUE if you only want to plot the significant hits from cpdb.
#' @param split.by column name in the metadata/coldata table to split the spots by. Can only take columns with binary options. If specified, name to split by MUST be specified in the meta file provided to cpdb prior to analysis.
#' @param gene.family default = NULL. some predefined group of genes. can take one of these options: "chemokines", "Th1", "Th2", "Th17", "Treg", "costimulatory", "coinhibitory", "niche"
#' @param genes default = NULL. can specify custom list of genes if gene.family is NULL
#' @param scale logical. scale the expression to mean +/- SD. NULL defaults to TRUE.
#' @param standard_scale logical. scale the expression to range from 0 to 1. NULL defaults to FALSE.
#' @param col_option specify plotting colours
#' @param noir default = FALSE. Ben's current phase. makes it b/w
#' @param highlight colour for highlighting p <0.05
#' @param ... passes arguments to grep for cell_type1 and cell_type2.
#' @return ggplot dot plot object of cellphone db output
#' @examples
#' scdata <- readRDS("./scdata.RDS", check.names = FALSE)
#' pvals <- read.delim("./cpdb/output/pvalues.txt", check.names = FALSE)
#' means <- read.delim("./cpdb/output/means.txt", check.names = FALSE)
#' plot_cpdb("Bcell", "Tcell", scdata, 'seurat_clusters', means, pvals, split.by = "group", genes = c("CXCL13", "CD274", "CXCR5"))
#' @import viridis
#' @import ggplot2
#' @import reshape2
#' @export

plot_cpdb <- function(cell_type1, cell_type2, scdata, idents, means_file, pvals_file, p.adjust.method = NULL, keep_significant_only = FALSE, split.by = NULL, gene.family = NULL, genes = NULL, scale = NULL, standard_scale = NULL, col_option = viridis::viridis(50), noir = FALSE, highlight = "red", ...) {
	require(ggplot2)
	require(reshape2)
	require(viridis)

	if (class(scdata) %in% c("SingleCellExperiment", "SummarizedExperiment")) {
		cat("data provided is a SingleCellExperiment/SummarizedExperiment object", sep = "\n")
		cat("extracting expression matrix", sep = "\n")
		require(SummarizedExperiment)
		require(SingleCellExperiment)
		exp_mat <- assay(scdata)
		metadata <- ColData(scdata)
	} else if (class(scdata) == "Seurat") {
		cat("data provided is a Seurat object", sep = "\n")
		cat("extracting expression matrix", sep = "\n")
		exp_mat <- tryCatch(scdata@data, error = function(e) {
			tryCatch(GetAssayData(object = scdata), error = function(e) {
				stop(paste0("are you sure that your data is normalized?"))
				return(NULL)
			})
		})
		metadata <- scdata@meta.data
	}

	means_mat <- means_file
	pvals_mat <- pvals_file
	rownames(means_mat) <- make.names(means_mat$interacting_pair, unique = TRUE)
	rownames(pvals_mat) <- make.names(pvals_mat$interacting_pair, unique = TRUE)
	colnames(means_mat) <- gsub("\\|", "-", colnames(means_mat))
	rownames(means_mat) <- gsub("_", "-", rownames(means_mat))
	rownames(means_mat) <- gsub("[.]", " ", rownames(means_mat))
	colnames(pvals_mat) <- gsub("\\|", "-", colnames(pvals_mat))
	rownames(pvals_mat) <- gsub("_", "-", rownames(pvals_mat))
	rownames(pvals_mat) <- gsub("[.]", " ", rownames(pvals_mat))

	if(length(p.adjust.method) > 0){
		pvals_tmp <- pvals[,12:ncol(pvals)]
		pvals_adj <- matrix(p.adjust(as.vector(as.matrix(pvals_tmp)), method=p.adjust.method),ncol=ncol(pvals_tmp))
		colnames(pvals_adj) <- colnames(pvals_tmp)
		pvals <- cbind(pvals[,c(1:11)], pvals_adj)
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

		labels <- unique(labels)
		groups <- unique(metadata[[split.by]])

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
			labels <- unique(labels)

			ct1 = grep(cell_type1, labels, value = TRUE, ...)
			ct2 = grep(cell_type2, labels, value = TRUE, ...)
		}
	} else {
		if(length(idents) > 1){
			labels <- idents
		} else {
			labels <- metadata[[idents]]
		}
		labels <- unique(labels)

		ct1 = grep(cell_type1, labels, value = TRUE, ...)
		ct2 = grep(cell_type2, labels, value = TRUE, ...)
		ct1 = paste0(ct1, collapse = '|')
		ct2 = paste0(ct2, collapse = '|')
	}

	if (ct1 == ''){ct1 = NA}
	if (ct2 == ''){ct2 = NA}
	checklabels2 <- any(colnames(means_mat) %in% c(ct1,ct2))

	if(!checklabels1){
		if(length(idents) > 1){
    		# relatively relaxed criteria to allow for the program to continue
			ct_1 <- grep(ct1, idents, value = TRUE, ...)
			ct_2 <- grep(ct2, idents, value = TRUE, ...)
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
			ct_1 <- grep(ct1, metadata[[idents]], value = TRUE, ...)
			ct_2 <- grep(ct2, metadata[[idents]], value = TRUE, ...)
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
		ct_1 <- grep(ct1, colnames(means_mat), value = TRUE)
		ct_2 <- grep(ct2, colnames(means_mat), value = TRUE)
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
		Th1 <- grep("IL2|IL12|IL18|IL27|IFNG|IL10|TNF$|TNF |LTA|LTB|STAT1|CCR5|CXCR3|IL12RB1|IFNGR1|TBX21|STAT4", means_mat$interacting_pair)
		Th2 <- grep("IL4|IL5|IL25|IL10|IL13|AREG|STAT6|GATA3|IL4R", means_mat$interacting_pair)
		Th17 <- grep("IL21|IL22|IL24|IL26|IL17A|IL17A|IL17F|IL17RA|IL10|RORC|RORA|STAT3|CCR4|CCR6|IL23RA|TGFB", means_mat$interacting_pair)
		Treg <- grep("IL35|IL10|FOXP3|IL2RA|TGFB", means_mat$interacting_pair)
		costimulatory <- grep("CD86|CD80|CD48|LILRB2|LILRB4|TNF|CD2|ICAM|SLAM|LT[AB]|NECTIN2|CD40|CD70|CD27|CD28|CD58|TSLP|PVR|CD44|CD55|CD[1-9]", means_mat$interacting_pair)
		coinhibitory <- grep("SIRP|CD47|ICOS|TIGIT|CTLA4|PDCD1|CD274|LAG3|HAVCR|VSIR", means_mat$interacting_pair)
		niche <- grep("CSF", means_mat$interacting_pair)
		query_group <- list(chemokines = chemokines, chemokine = chemokines, Th1 = Th1, Th2 = Th2, Th17 = Th17, Treg = Treg, costimulatory = costimulatory, coinhibitory = coinhibitory, costimulation = costimulatory, coinhibition = coinhibitory, niche = niche)
	} else if (is.null(gene.family) & !is.null(genes)){
		query <- grep(paste(genes, collapse="|"), means_mat$interacting_pair)
	}

	create_celltype_query <- function(ctype1, ctype2){
		ct1 = list()
		ct2 = list()
		for (i in 1:length(ctype2)){
			ct1[i] = paste0(ctype1, "-", ctype2[i])
			ct2[i] = paste0(ctype2[i], "-", ctype1)
		}
		ct_1 = do.call(paste0, list(ct1, collapse = '|'))
		ct_2 = do.call(paste0, list(ct2, collapse = '|'))
		ct = list(ct_1, ct_2)
		ct = do.call(paste0, list(ct, collapse = '|'))
		return(ct)
	}

	keep_interested_groups <- function(g, ct){
		ctx <- strsplit(ct, "\\|")[[1]]
		idx <- grep(paste0(g, ".*-", g), ctx)
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

		labels <- unique(labels)
		groups <- unique(metadata[[split.by]])

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
				celltype[[i]] <- create_celltype_query(c_type1[[i]], c_type2)
				celltype[[i]] <- lapply(grp, keep_interested_groups, celltype[[i]])
			}

			for (i in 1:length(celltype)){
				celltype[[i]] <- celltype[[i]][-which(celltype[[i]] == '')]
			}

			celltype <- lapply(celltype, unlist)
			cell_type <- do.call(paste0, list(celltype, collapse = "|"))
		} else {
			if(length(idents) > 1){
				labels <- idents
			} else {
				labels <- metadata[[idents]]
			}
			labels <- unique(labels)

			c_type1 = as.list(grep(cell_type1, labels, value = TRUE, ...))
			c_type2 = as.list(grep(cell_type2, labels, value = TRUE, ...))

			celltype = list()
			for (i in 1:length(c_type1)){
				celltype[[i]] <- create_celltype_query(c_type1[[i]], c_type2)
			}
			cell_type <- do.call(paste0, list(celltype, collapse = "|"))
		}
	} else {
		if(length(idents) > 1){
			labels <- idents
		} else {
			labels <- metadata[[idents]]
		}
		labels <- unique(labels)

		c_type1 = as.list(grep(cell_type1, labels, value = TRUE, ...))
		c_type2 = as.list(grep(cell_type2, labels, value = TRUE, ...))

		celltype = list()
		for (i in 1:length(c_type1)){
			celltype[[i]] <- create_celltype_query(c_type1[[i]], c_type2)
		}
		cell_type <- do.call(paste0, list(celltype, collapse = "|"))
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
			pvals_mat <- pvals_mat[,as.vector(group_id)]
		}
	}

	if(nrow(means_mat) > 2){
		d <- dist(as.data.frame(means_mat))
		h <- hclust(d)
		means_mat <- means_mat[h$order, ]
		pvals_mat <- pvals_mat[h$order, ]
	}

	# scaling
	range01 <- function(x){(x-min(x))/(max(x)-min(x))}
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
	means_mat2[which(means_mat == 0)] <- NA

	# remove rows that are entirely NA
	pvals_mat2 <- pvals_mat2[rowSums(is.na(means_mat2)) != ncol(means_mat2), ,drop = FALSE]
	means_mat2 <- means_mat2[rowSums(is.na(means_mat2)) != ncol(means_mat2), ,drop = FALSE]

	if ((length(standard_scale) > 0 && standard_scale) | (length(scale) > 0 && scale) | (length(scale) < 1 && length(standard_scale) < 1)){
		df_means <- melt(means_mat2, value.name = "scaled_means")
	} else {
		df_means <- melt(means_mat2, value.name = "means")
	}

	if(length(p.adjust.method) > 0){
		df_pvals <- melt(pvals_mat2, value.name = "padj")
		df <- data.frame(cbind(df_means, padj = df_pvals$padj))
		df$padj[which(df$padj >= 0.05)] <- NA
		if (keep_significant_only){
			# keep the entire row/ all the comparisons
			df_ <- split(df, as.character(df$Var1))
			anysig <- lapply(df_, function(x){
				keep <- any(x$padj < 0.05)
				return(keep)
			})
			df_ <- df_[which(unlist(anysig))]
			names(df_) <- NULL
			df <- do.call(rbind, df_)			
		}
		df$padj[which(df$padj == 0)] <- 0.001
	} else {
		df_pvals <- melt(pvals_mat2, value.name = "pvals")
		df <- data.frame(cbind(df_means, pvals = df_pvals$pvals))
		df$pvals[which(df$pvals >= 0.05)] <- NA
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
	
	if ((length(standard_scale) > 0 && standard_scale) | (length(scale) > 0 && scale) | (length(scale) < 1 && length(standard_scale) < 1)){
		if(length(p.adjust.method) > 0 && p.adjust.method != 'none'){
			g <- ggplot(df, aes(x = Var2, y = Var1, color = -log10(padj), fill = scaled_means, size = scaled_means))
		} else {
			g <- ggplot(df, aes(x = Var2, y = Var1, color = -log10(pvals), fill = scaled_means, size = scaled_means))
		}
	} else {
		if(length(p.adjust.method) > 0 && p.adjust.method != 'none'){
			g <- ggplot(df, aes(x = Var2, y = Var1, color = -log10(padj), fill = means, size = means))
		} else {
			g <- ggplot(df, aes(x = Var2, y = Var1, color = -log10(pvals), fill = means, size = means))
		}
	}

	g <- g + geom_point(pch = 21) +
	theme_bw() +
	theme(axis.text.x = element_text(angle = 45, hjust = 0),
		axis.ticks = element_blank(),
		axis.title.x = element_blank(),
		axis.title.y = element_blank()) +
	scale_x_discrete(position = "top") +
	scale_color_gradientn(colors = highlight, na.value = "white") +
	scale_radius(range = c(0,5))

	if(noir){
		g <- g + scale_fill_gradient(low = "white", high = "#131313", na.value = "white")
	} else {
		if(length(col_option) == 1){
			g <- g + scale_fill_gradientn(colors = colorRampPalette(c("white", col_option))(100), na.value = "white")
		} else {
			g <- g + scale_fill_gradientn(colors = c("white", colorRampPalette(col_option)(99)), na.value = "white")
		}
	}
	if(!is.null(gene.family) & is.null(genes)){
		g <- g + ggtitle(gene.family)
	}
	return(g)
}
