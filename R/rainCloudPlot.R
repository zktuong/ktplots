#' rainCloudPlot
#'
#' @param data plotting data
#' @param grouby x axis
#' @param parameter y axis
#' @param default default plotting options
#' @param ... passed to geom_glat_violin
#' @return rainCloudPlot
#' @examples
#' rainCloudPlot(diamonds, cut, carat) + coord_flip()
#' @import dplyr
#' @import ggplot2 
#' @export
#'
rainCloudPlot <- function(data, groupby, parameter, default = TRUE, ...){

	g <- ggplot(data, aes(x = get(grouby), y = get(parameter), fill = get(groupby))) +

	if (default){
		g <- g + geom_flat_violin(position = position_nudge(x = .1, y = 0), adjust = 1.5, trim = FALSE, alpha = .5, colour = NA, ...) +	
		geom_point(aes(x = as.numeric(get(groupby))-.15, y = get(parameter), colour = get(groupby)),position = position_jitter(width = .05), size = 1, shape = 20) +
		geom_boxplot(outlier.shape = NA, alpha = .5, width = .1, colour = "black")
	} else {
		g <- g + geom_flat_violin(...)+
		geom_point(aes(x = as.numeric(get(groupby))-.15, y = get(parameter), colour = get(groupby)),position = position_jitter(width = .05), size = 1, shape = 20) +
		geom_boxplot(outlier.shape = NA, alpha = .5, width = .1, colour = "black")
	}

	return(g)
}