#' Plotting cellphonedb results
#' 
#' @param cell_type1 cell type 1
#' @param cell_type2 cell type 2
#' @param means_file object holding means.txt from cpdb output
#' @param pvals_file object holding pvals.txt from cpdb output
#' @param groups conditions to do the plotting
#' @param gene.family default = NULL. some predefined group of genes. can take one of these options: "chemokines", "Th1", "Th2", "Th17", "Treg", "costimulatory", "coinhibitory", "niche"
#' @param genes default = NULL. can specify custom list of genes if gene.family is NULL
#' @return ggplot dot plot object of cellphone db output
#' @examples
#' pvals <- read.delim("./cpdb/output/pvalues.txt", check.names = FALSE)
#' means <- read.delim("./cpdb/output/means.txt", check.names = FALSE) 
#' plot_cpdb("MNP_5", "Lymphoid_[12345]", means, pvals, groups = c("normal", "tumor"), genes = c("CXCL13", "CD274", "CXCR5"))
#' @import viridis
#' @import ggplot2
#' @import reshape2
#' @export

plot_cpdb <- function(cell_type1, cell_type2, means_file, pvals_file, groups, gene.family = NULL, genes = NULL, ...) {
	require(ggplot2)
	require(reshape2)
	require(viridis)

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
	
	if (is.null(gene.family) & is.null(genes)){
		stop("Please specify either genes or gene.family")
		cat("gene.family can be one of the following:", sep = "\n")
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
		query_group <- list(chemokines = chemokines, Th1 = Th1, Th2 = Th2, Th17 = Th17, Treg = Treg, costimulatory = costimulatory, coinhibitory = coinhibitory, niche = niche)
	} else if (is.null(gene.family) & !is.null(genes)){
		query <- grep(paste(genes, collapse="|"), means_mat$interacting_pair)
	} 

	group1 = groups[1]
	group2 = groups[2]

	cell_type = paste0(paste0(group1, ".*", cell_type1, ".*-", group1, ".*", cell_type2), "|", paste0(group1,".*", cell_type2, ".*-", group1,".*", cell_type1), "|", paste0(group2,".*", cell_type1, ".*-", group2, ".*", cell_type2), "|", paste0(group2,".*", cell_type2, ".*-", group2, ".*", cell_type1)) 

	if(!is.null(gene.family) & is.null(genes)){
		means_mat <- means_mat[query_group[[gene.family]], grep(cell_type, colnames(means_mat))]
		pvals_mat <- pvals_mat[query_group[[gene.family]], grep(cell_type, colnames(pvals_mat))]
	} else if (is.null(gene.family) & !is.null(genes)){
		means_mat <- means_mat[query, grep(cell_type, colnames(means_mat))]
		pvals_mat <- pvals_mat[query, grep(cell_type, colnames(pvals_mat))]
	}
	
	# rearrange the columns so that it interleaves the two groups
	group.1 <- grep(group1, colnames(means_mat))
	group.2 <- grep(group2, colnames(means_mat))

	means_mat <- means_mat[,as.vector(rbind(group.1, group.2))]
	pvals_mat <- pvals_mat[,as.vector(rbind(group.1, group.2))]

	if(nrow(means_mat) > 2){
		d <- dist(as.data.frame(means_mat))
		h <- hclust(d)
		means_mat <- means_mat[h$order, ]
		pvals_mat <- pvals_mat[h$order, ]		
	} 
	
	# scale
	means_mat2 <- t(scale(t(means_mat)))
	pvals_mat2 <- as.matrix(pvals_mat)
	means_mat2 <- as.matrix(means_mat2)
	means_mat2[which(means_mat == 0)] <- NA

	# remove rows that are entirely NA
	pvals_mat2 <- pvals_mat2[rowSums(is.na(means_mat2)) != ncol(means_mat2), ,drop = FALSE]
	means_mat2 <- means_mat2[rowSums(is.na(means_mat2)) != ncol(means_mat2), ,drop = FALSE]

	df_means <- melt(means_mat2, value.name = "scaled_means")
	df_pvals <- melt(pvals_mat2, value.name = "pvals")
	df <- data.frame(cbind(df_means, pvals = df_pvals$pvals))
	df$pvals[which(df$pvals >= 0.05)] <- NA
	df$pvals[which(df$pvals == 0)] <- 0.001
	df$group <- gsub(paste0(group1,"_|", group2, "_"), "", df$Var2)

	g <- ggplot(df, aes(x = Var2, y = Var1, color = -log10(pvals), fill = scaled_means, size = scaled_means)) + 
  		geom_point(pch = 21) +
  		theme_bw() +
  		theme(axis.text.x = element_text(angle = 45, hjust = 0),
  			axis.ticks = element_blank(),
  			axis.title.x = element_blank(),
  			axis.title.y = element_blank()) +
  		scale_x_discrete(position = "top") +
  		scale_color_gradientn(colors = "red", na.value = "white") +
  		scale_fill_gradientn(colors = c("white", viridis(50, direction = 1)), na.value = "white") +
  		scale_size_continuous(range = c(0,5))
  	if(!is.null(gene.family) & is.null(genes)){
  		g <- g + ggtitle(gene.family)
  	}
  	return(g)
}