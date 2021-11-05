#' Plotting cellphonedb results
#'
#' @param cell_type1 cell type 1
#' @param cell_type2 cell type 2
#' @param scdata single-cell data. can be seurat/summarizedexperiment object
#' @param idents vector holding the idents for each cell or column name of scdata's metadata. MUST match cpdb's columns
#' @param means object holding means.txt from cpdb output
#' @param pvals object holding pvals.txt from cpdb output
#' @param deconvoluted object holding deconvoluted.txt from cpdb output
#' @param p.adjust.method correction method. p.adjust.methods of one of these options: c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none")
#' @param keep_significant_only logical. Default is FALSE. Switch to TRUE if you only want to plot the significant hits from cpdb.
#' @param split.by column name in the metadata/coldata table to split the spots by. Can only take columns with binary options. If specified, name to split by MUST be specified in the meta file provided to cpdb prior to analysis.
#' @param scale logical. scale the expression to mean +/- SD. NULL defaults to TRUE.
#' @param standard_scale logical. scale the expression to range from 0 to 1. Default is TRUE
#' @param separator default = NULL. separator to use to split between celltypes. Unless otherwise specified, the separator will be `>@<`. Make sure the idents and split.by doesn't overlap with this.
#' @param gene_symbol_mapping default = NULL.column name for rowData in sce holding the actual gene symbols if row names aren't gene symbols
#' @param frac default = 0.2. Does nothing at the moment.
#' @param remove_self default = TRUE. Remove self-self arcs. might not be working.
#' @param desiredInteractions default = NULL. Specific list of celltype comparisons e.g. list(c('CD4_Tcm', 'cDC1'), c('CD4_Tcm', 'cDC2')). Also accepts a dataframe where first column is celltype 1 and 2nd column is celltype 2
#' @param interaction_grouping default = NULL. dataframe specifying groupings of cellphonedb interactions. First column must be cellphonedb's interacting_pair column. second column is whatever grouping you want.
#' @param edge_group_colors default = NULL. vector for colour mapping for edge groups. only used if split.by is specified.
#' @param node_group_colors default = NULL. vector for colour mapping for node labels.
#' @param ... passes arguments plot_cpdb
#' @return ggplot dot plot object of cellphone db output
#' @examples
#' \donttest{
#' }
#' @import ggplot2
#' @import ggraph
#' @import ggrepel
#' @export

plot_cpdb2 <- function(cell_type1, cell_type2, scdata, idents, means, pvals, deconvoluted, p.adjust.method = NULL, keep_significant_only = TRUE, split.by = NULL, standard_scale = TRUE, separator = NULL, gene_symbol_mapping = NULL, frac = 0.2, remove_self = TRUE, desiredInteractions = NULL, interaction_grouping = NULL, edge_group_colors = NULL, node_group_colors = NULL, ...){

	if (length(separator) > 0){
		sep = separator
	} else {
		sep = '>@<'
	}	

	cpdb_int = plot_cpdb(cell_type1 = cell_type1, cell_type2 = cell_type2, scdata = scdata, idents = idents, split.by = split.by, means = means, pvals = pvals, keep_significant_only = keep_significant_only, standard_scale = standard_scale, ...)
	requireNamespace('SummarizedExperiment')
	requireNamespace('SingleCellExperiment')
	lr_interactions <- cpdb_int$data
	subset_clusters <- unique(unlist(lapply(as.list(lr_interactions$group), strsplit, sep)))
	sce_subset <- scdata[, SummarizedExperiment::colData(scdata)[,idents] %in% subset_clusters]
	interactions <- means[,c('interacting_pair', 'gene_a', 'gene_b', 'partner_a', 'partner_b')]
	interactions$converted <- gsub('-', ' ', interactions$interacting_pair)
	interactions$converted <- gsub('_', '-', interactions$converted)
	interactions_subset <- interactions[interactions$converted %in% lr_interactions$Var1, ]
	if (any(grepl('ENSG', interactions_subset$gene_a))){
		tm0 <- do.call(c, lapply(as.list(interactions_subset$interacting_pair), strsplit, '_'))
		tm0 = data.frame(t(matrix(unlist(tm0),2, length(unlist(tm0))/2)))
		colnames(tm0) <- c('id_a', 'id_b')
		interactions_subset <- cbind(interactions_subset, tm0)
		dictionary_a <- interactions_subset$id_a
		names(dictionary_a) <- interactions_subset$gene_a
		dictionary_b <- interactions_subset$id_b
		names(dictionary_b) <- interactions_subset$gene_b
		dictionary <- c(as.list(dictionary_a), as.list(dictionary_b))
		if (length(which(names(dictionary) == '')) > 0){
			dictionary <- dictionary[-which(names(dictionary) == '')]
		}
	}
	
	if (!is.null(interaction_grouping)){
		if ((class(interaction_grouping) == 'data.frame')){
		interactions_subset$group = interaction_grouping[,2][match(interactions_subset$interacting_pair, 	interaction_grouping[,1])]
		}
	}
	
	# extract all the possible genes.
	geneid = unique(c(interactions_subset$gene_a, interactions_subset$gene_b))
	rmg = which(geneid == "")
	if (length(rmg) > 0){
		geneid = geneid[-which(geneid == "")]	
	}
	sce_subset_tmp <- sce_subset[geneid, ]
	
	# split to list and calculate celltype mean for each treatment group
	if (class(scdata) %in% c("SingleCellExperiment", "SummarizedExperiment")) {		
		requireNamespace('SummarizedExperiment')
		requireNamespace('SingleCellExperiment')
		meta <- as.data.frame(SummarizedExperiment::colData(sce_subset_tmp))
	} else if (class(scdata) == "Seurat") {
		stop("sorry not supported yet")
	}

	sce_list <- list()
	sce_list_alt <- list()
	for (x in unique(meta[,split.by])){
		sce_list[[x]] <- list()
		sce_list_alt[[x]] <- list()
	}

	for (n in names(sce_list)){
		for (x in unique(meta[,idents])){
			sce_list[[n]][[x]] <- sce_subset_tmp[, meta[,idents] == x & meta[,split.by] == n]
			sce_list_alt[[n]][[x]] <- sce_subset[, meta[,idents] == x & meta[,split.by] == n]
		}
	}
	
	cellTypeMeans <- function(x){
		cm <- Matrix::rowMeans(SingleCellExperiment::counts(x))
		return(cm)
	}
	
	cellTypeFraction <- function(x){
		cm <- Matrix::rowMeans(SingleCellExperiment::counts(x) > 0)
		return(cm)
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
	
	for (n in names(sce_list2)){
		colnames(sce_list2[[n]]) <- paste0(paste0(n, '_') , colnames(sce_list2[[n]]))
		colnames(sce_list3[[n]]) <- paste0(paste0(n, '_') , colnames(sce_list3[[n]]))
	}

	sce_list2 <- do.call(cbind, sce_list2)
	sce_list3 <- do.call(cbind, sce_list3)
	
	humanreadablename = c()
	for (i in row.names(sce_list2)){
		humanreadablename = c(humanreadablename, as.character(unlist(dictionary[i])))
	}
	
	rownames(sce_list2) <- humanreadablename
	rownames(sce_list3) <- humanreadablename
	
	findComplex <- function(interaction){
		idxa <- which(interaction$gene_a == '')
		idxb <- which(interaction$gene_b == '')
		complexa <- gsub('complex:', '', interaction$partner_a[idxa])
		complexb <- gsub('complex:', '', interaction$partner_b[idxb])
	
		if (length(complexa) > 0){
			if (length(complexb) > 0) {
				res <- c(complexa, complexb)
			} else {
				res <- complexa
			}
		} else if (length(complexb) > 0){
			res <- complexb
		} else {
			res <- NULL
		}
		return(res)
	}
	
	decon_subset <- deconvoluted[deconvoluted$complex_name %in% findComplex(interactions_subset), ]
	if (nrow(decon_subset) > 0){
		# although multiple rows are returned, really it's the same value for the same complex
		decon_subset <- split(decon_subset, decon_subset$complex_name)
		decon_subset_expr <- lapply(decon_subset, function(x) {
			x <- x[,colnames(sce_list2)]
			x <- colMeans(x)
			return(x)
		})

		# to retrieve the fraction, we use the average fraction of all genes mapping to the complex
		cellTypeFraction_complex <- function(sce_, genes, gene_symbol_mapping = NULL){
			scex <- tryCatch(sce_[genes, ], error = function(e){
				if (!is.null(gene_symbol_mapping)){
					sce_[which(rowData(sce_)[,gene_symbol_mapping] %in% genes), ]
				} else {
					sce_[which(rowData(sce_)[,'index'] %in% genes), ]
				}
			})
		
			cm <- mean(Matrix::rowMeans(SingleCellExperiment::counts(scex) > 0))
			return(cm)
		}
		
		decon_subset_fraction <- lapply(decon_subset, function(x) {
			x <- unique(x$gene_name)
			test <- lapply(sce_list_alt, function(y){
				return(lapply(y, cellTypeFraction_complex, x, gene_symbol_mapping))
			})
			return(test)
		})
		
		decon_subset_fraction <- lapply(decon_subset_fraction, function(x) {
			y <- lapply(x, function(z) do.call(cbind, z))
			for (i in 1:length(y)){
				colnames(y[[i]]) <- paste0(names(y[i]), '_', colnames(y[[i]]))
			}
			y <- do.call(cbind, y)
			return(y)
		})
		
		decon_subset_expr <- do.call(rbind, decon_subset_expr)
		decon_subset_fraction <- do.call(rbind, decon_subset_fraction)
		row.names(decon_subset_fraction) <- row.names(decon_subset_expr)
		
		expr_df <- rbind(sce_list2, decon_subset_expr)
		fraction_df <- rbind(sce_list3, decon_subset_fraction)
	} else {
		expr_df <- sce_list2
		fraction_df <- sce_list3
	}
	

	# make a big fat edgelist
	if (!is.null(desiredInteractions)){
		if (class(desiredInteractions) == 'list'){
			desiredInteractions_ <- c(desiredInteractions, lapply(desiredInteractions, rev))
			cell_type_grid <- as.data.frame(do.call(rbind, desiredInteractions_))
		} else if ((class(desiredInteractions) == 'data.frame')){
			cell_type_grid <- desiredInteractions
		}
		cells_test = unique(unlist(desiredInteractions))
	} else {
		cells_test <- unique(droplevels(meta[,idents]))
		cell_type_grid <- expand.grid(cells_test, cells_test)
	}
	
	if (remove_self){
		rm_idx <- which(cell_type_grid[,1] == cell_type_grid[,2])
		if (length(rm_idx) > 0){
			cell_type_grid <- cell_type_grid[-rm_idx, ]
		}
	}
	
	ligand <- interactions_subset$id_a
	receptor <- interactions_subset$id_b
	pair <- interactions_subset$interacting_pair
	converted_pair <- interactions_subset$converted
	producers <- as.character(cell_type_grid[, 1])
	receivers <- as.character(cell_type_grid[, 2])

	generateDf <- function(ligand, receptor, pair, converted_pair, producers, receivers, cell_type_means, cell_type_fractions, splitted = NULL){
	if (!is.null(splitted)){
		pp <- paste0(splitted, '_', producers)
		rc <- paste0(splitted, '_', receivers)
	} else {
		pp <- producers
		rc <- receivers
	}
	df_ <- data.frame(
		'ligand' = ligand,
		'receptor' = receptor,
		'pair' = pair,
		'producer' = pp,
		'receiver' = rc,
		'producer_expression' = as.numeric(cell_type_means[ligand , pp]),
		'producer_fraction' = as.numeric(cell_type_fractions[ligand, pp]),
		'receiver_expression' = as.numeric(cell_type_means[receptor, rc]),
		'receiver_fraction' = as.numeric(cell_type_fractions[receptor,rc]))
	df_$from = paste0(df_$producer, "_", df_$ligand)
	df_$to = paste0(df_$receiver, "_", df_$receptor)
	if (!is.null(splitted)){		
		df_$producer_ = df_$producer
		df_$receiver_ = df_$receiver
		df_$from = gsub(paste0(splitted, '_'), '', df_$from)
		df_$to = gsub(paste0(splitted, '_'), '', df_$to)
		df_$producer = gsub(paste0(splitted, '_'), '', df_$producer)
		df_$receiver = gsub(paste0(splitted, '_'), '', df_$receiver)
		df_$barcode = paste0(df_$producer_, '-', df_$receiver_, sep, converted_pair)
	} else {
		df_$barcode = paste0(df_$producer, '-', df_$receiver, sep, converted_pair)
	}

	return(df_)}

	dfx <- list()
	if (!is.null(split.by)){
		for (i in unique(meta[,split.by])){
			dfx[[i]] = generateDf(ligand, receptor, pair, converted_pair, producers, receivers, expr_df, fraction_df, i)
		}
	} else {
		dfx[[1]] = generateDf(ligand, receptor, pair, converted_pair, producers, receivers, expr_df, fraction_df)
	}

#set the bundled connections
	df0 <- lapply(dfx, function(x) x[x$producer_fraction > frac & x$receiver_fraction > frac, ]) #save this for later

# now construct the hierachy
	constructGraph <- function(el, el0, unique_id, interactions_df, plot_cpdb_out, edge_group = FALSE, edge_group_colors = NULL,  node_group_colors = NULL){
		require(igraph)
		celltypes <- unique(c(as.character(el$producer), as.character(el$receiver)))
		el1 <- data.frame(from = "root", to = celltypes, barcode_1 = NA, barcode_2 = NA, barcode_3 = NA)
		el2 <- data.frame(from = celltypes, to = paste0(celltypes, "_", "ligand"), barcode_1 = NA, barcode_2 = NA, barcode_3 = NA)
		el3 <- data.frame(from = celltypes, to = paste0(celltypes, "_", "receptor"), barcode_1 = NA, barcode_2 = NA, barcode_3 = NA)
		el4 <- do.call(rbind, lapply(celltypes, function(x){
			cell_ligands <- grep(x, el$from, value = TRUE)
			cell_ligands_idx <- grep(x, el$from)
			if(length(cell_ligands) >0){
				df <-   data.frame(from = paste0(x, "_", "ligand"), to = cell_ligands, barcode_1 = el$barcode[cell_ligands_idx], barcode_2 = el$pair[cell_ligands_idx], barcode_3 = paste0(el$from[cell_ligands_idx], sep, el$to[cell_ligands_idx]))
			} else { df = NULL}}))
		el5 <- do.call(rbind, lapply(celltypes, function(x){
			cell_ligands <- grep(x, el$to, value = TRUE)
			cell_ligands_idx <- grep(x, el$to)
			if(length(cell_ligands) >0){
				df <-   data.frame(from = paste0(x, "_", "receptor"), to = cell_ligands, barcode_1 = el$barcode[cell_ligands_idx], barcode_2 = el$pair[cell_ligands_idx], barcode_3 = paste0(el$from[cell_ligands_idx], sep, el$to[cell_ligands_idx]))
			}	else{df = NULL}
		}))

		gr_el <- do.call(rbind, list(el1, el2, el3, el4, el5))
		plot_cpdb_out$barcode <- paste0(plot_cpdb_out$Var2, sep, plot_cpdb_out$Var1)
		mean_col <- grep('means$', colnames(plot_cpdb_out), value = TRUE)
		means <- plot_cpdb_out[match(gr_el$barcode_1, plot_cpdb_out$barcode), mean_col]
		pval_col <- grep('pvals', colnames(plot_cpdb_out), value = TRUE)
		pvals <- plot_cpdb_out[match(gr_el$barcode_1, plot_cpdb_out$barcode), pval_col]
		gr_el <- cbind(gr_el, means, pvals)

		if (edge_group){
			groups <- interactions_df$group[match(gr_el$barcode_2, interactions_df$interacting_pair)]
		}

		gr <- graph_from_edgelist(as.matrix(gr_el[,1:2]))
		E(gr)$interaction_score <- means
		E(gr)$pvals <- pvals
		if (edge_group){
			E(gr)$group <- groups
		}
		E(gr)$name <- gr_el$barcode_3

	# order the graph vertices
		V(gr)$type <- NA
		V(gr)$type[V(gr)$name %in%  el4$to] <- "ligand"
		V(gr)$type[V(gr)$name %in%  el5$to] <- "receptor"

		from = match(el0$from, V(gr)$name)
		to = match(el0$to, V(gr)$name)
		dat = data.frame(from = el0$from, to = el0$to)
		dat$barcode = paste0(dat$from, sep, dat$to)
		interaction_score = E(gr)$interaction_score[match(dat$barcode, gr_el$barcode_3)]
		pval = E(gr)$pvals[match(dat$barcode, gr_el$barcode_3)]
		pval[is.na(pval)] <- 1
		pval <- range01(-log10(pval))

		if (edge_group){
			group = E(gr)$group[match(dat$barcode, gr_el$barcode_3)]
		}

		ligand_expr <- data.frame("cell_mol" = el$from,
			"expression" = el$producer_expression,
			"fraction" = el$producer_fraction)
		recep_expr <- data.frame("cell_mol" = el$to,
			"expression" = el$receiver_expression,
			"fraction" = el$receiver_fraction)
		expression <- rbind(ligand_expr, recep_expr)

		V(gr)$expression = NA
		V(gr)$expression[match(expression$cell_mol, V(gr)$name)] <- expression$expression

		V(gr)$fraction = 0
		V(gr)$fraction[match(expression$cell_mol, V(gr)$name)] <- expression$fraction

  	#celltype
		V(gr)$celltype <- NA
		for(x in cells_test){
			idx <- grepl(x, V(gr)$name)
			V(gr)$celltype[idx] <- x
		}

	# V(gr)$name[V(gr)$fraction < frac] = NA
		V(gr)$name[!V(gr)$name %in% c(el0$from, el0$to)] <- NA
	# final adjustments to the name

		for(x in unique_id){
			V(gr)$name <- gsub(paste0(x, '_'), '', V(gr)$name)
		}

		library(ggraph)
		library(ggrepel)

		if (!is.null(edge_group_colors)){
			edge_group_colors = edge_group_colors
		} else {
			nn = length(unique(E(gr)$group))
			edge_group_colors = kelvinny::gg_color_hue(nn)
		}
		if (!is.null(node_group_colors)){
			node_group_colors = node_group_colors
		} else {
			nn = length(unique(meta[,idents]))
			node_group_colors = kelvinny::gg_color_hue(nn)
		}
  	#plot the graph
		pl <- ggraph(gr, layout = 'dendrogram', circular = TRUE) +
		geom_conn_bundle(data = get_con(from = from, to = to, group = group, `-log10(sig)` = pval, interaction_score = interaction_score),
			aes(colour = group, alpha = interaction_score, width = `-log10(sig)`), tension = 0.5) +
		scale_edge_width(range = c(1, 3)) +
		scale_edge_alpha(range = c(0, 3), limits = c(0, 1)) +
		scale_edge_color_manual(values = edge_group_colors)+
		geom_node_point(pch =19, aes(size = fraction, filter = leaf, color = celltype, alpha = type)) +
		theme_void() + coord_fixed() +
		scale_size_continuous(limits= c(0, 1))  +
		scale_shape_manual(values = c("ligand" = 19, "receptor" = 15)) +
		scale_color_manual(values =  node_group_colors) +
		geom_text_repel(aes(x=x, y=y, label = name), segment.square = TRUE, segment.inflect = TRUE, segment.size = 0.2, force=0.5, size = 2, force_pull = 0) +
    	# geom_node_text(aes(x = x*1.15, y=y*1.15, filter = leaf, label=name, size =0.1))  +
		scale_alpha_manual(values = c("ligand" = 0.5,
			"receptor" = 1)) +
		small_legend(keysize = 0.5)
		return(pl)
	}

	gl <- list()
	if (!is.null(split.by)){
		edge_group = TRUE
	} else {
		edge_group = FALSE
	}

	for (i in 1:length(dfx)){
		gl[[i]] <- constructGraph(dfx[[i]], df0[[i]], cells_test, interactions_subset, lr_interactions, edge_group, edge_group_colors, node_group_colors)
	}

	if (length(gl) > 1){	
		return(gl)
	} else {
		return(gl[[1]])	
	}
}


