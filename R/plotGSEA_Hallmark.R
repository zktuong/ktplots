#' Plotting GSEA enrichment as dotplots
#' 
#' @param gsea processed gsea data table
#' @param group_ref if more than 2 comparisons, specify the name of the group to act as reference. Otherwise, use 1 or 2 to specify the order if the reference. Defaults to using the second level.
#' @param cols replacement colours
#' @param newlabels replacement labels
#' @return ggplot dot plot object of functions
#' @examples
#' plotGSEA_Hallmark(gsea_result) + ggtitle("GSEA_Hallmarks")
#' @import dplyr
#' @import Matrix
#' @import ggplot2
#' @import reshape2
#' @export

plotGSEA_Hallmark <- function(gsea, group_ref = NULL, cols = NULL, newlabels = NULL) {
	require(ggplot2)
	gsea$NES[which(is.na(gsea$NES))] <- 0
	gsea$ranking[which(is.na(gsea$ranking))] <- 0
	gsea <- gsea[order(gsea$ranking),]		
	gsea_spl <- split(gsea, gsea$group)
	if(!is.null(group_ref)){
		gsea_spl[[group_ref]] <- gsea_spl[[group_ref]][order(gsea_spl[[group_ref]]$ranking),]
		gsea_spl[[group_ref]]$ranking <- gsea_spl[[group_ref]]$ranking*999
	} else {
		gsea_spl[[2]] <- gsea_spl[[2]]$ranking*999
	}
	names(gsea_spl) <- NULL

	gsea <- do.call(rbind, gsea_spl)
	gsea <- gsea[order(gsea$ranking), ]
	gsea$pathway <- gsub("HALLMARK_|", "", gsea$pathway)
	gsea$group[which(gsea$pval >= 0.05 & gsea$padj >= 0.25)] <- "NA"
	
	x_lim_min <- abs(ceiling(min(-log10(gsea$pval))))
	x_lim_max <- abs(ceiling(max(-log10(gsea$pval))))

	if(x_lim_min > x_lim_max){
		xval1 <- x_lim_min * -1
		xval2 <- x_lim_min
	} else {
		xval1 <- x_lim_max * -1
		xval2 <- x_lim_max
	}

	if(!is.null(cols)){
		gg_color_hue <- function(n) {
			hues = seq(15, 375, length = n + 1)
			hcl(h = hues, l = 65, c = 100)[1:n]
		}
		cols. = gg_color_hue(length(dplyr::n_distict(gsea$group, na.rm = TRUE)))
	} else {
		cols. = cols
	}

	g <- ggplot(gsea, aes(x = -log10(pval)*sign(NES), y = reorder(pathway, ranking), col = group, size = NES)) + 
		geom_point() + 
		labs(x = expression(paste("Signed", " -log" ["10"], "P-value")), y = "Hallmarks") +
		theme_bw() +
		geom_vline(xintercept = 0) +
		geom_vline(xintercept = -log10(0.25)) +
		geom_vline(xintercept = -log10(0.25)*-1) +
		xlim(xval1, xval2) +
		scale_size_continuous(range = c(1,4)) +
		theme(panel.grid.major = element_blank(), 
			panel.grid.minor = element_blank(), 
			panel.background = element_blank(), 
			axis.line = element_blank(), 
			axis.ticks = element_blank())
	if(!is.null(newlabels))
		g <- g + scale_color_manual(values = cols., na.value = 'grey90', drop = FALSE, labels = newlabels)
	else {
		g <- g + scale_color_manual(values = cols., na.value = 'grey90', drop = FALSE)
	}
	g$data <- g$data[order(g$data$group, na.last = FALSE), ]
	return(g)
}

