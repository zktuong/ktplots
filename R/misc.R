#' miscellaneous functions
#' 
#' @return miscellaneous functions
#' @import dplyr
#' @import ggplot2 

#' @rdname nin
#' @name ktplots
#' @usage a \%nin\% b
#' @param a a valute
#' @param b a value 
#' @export
"%nin%" <- function(a, b) {
  if (!is.null(a)) a else b
}

#' @name ktplots
#' @param fontsize float/int
#' @param ... passed to ggplot2::theme
#' @export
small_legend <- function(fontsize = 5, ...){
	small_legend_theme <- theme(
		legend.title = element_text(size = fontsize), 
		legend.text  = element_text(size = fontsize),
		legend.key.size = unit(0.1, "lines"), ...)
	return(small_legend_theme)
}

#' @name ktplots
#' @param guidesize float/int
#' @param ... passed to ggplot2::theme
#' @export
small_guide <- function(guidesize = 1, ...){
	small_guide <- guides( 
		shape = guide_legend(override.aes = list(size = guidesize)), 
		color = guide_legend(override.aes = list(size = guidesize)), ...)
	return(small_guide)
}

#' @name ktplots
#' @param legendmargin margin(float/int, float/int, float/int, float/int)
#' @param ... passed to ggplot2::theme
#' @export
topright_legend <- function(legendmargin = margin(6, 6, 6, 6), ...){
	legend <- theme(legend.position = c(.99, .99), legend.justification = c('right', 'top'), legend.box.just = "right", legend.margin = legendmargin, ...)
	return(legend)
}

#' @name ktplots
#' @param legendmargin margin(float/int, float/int, float/int, float/int)
#' @param ... passed to ggplot2::theme
#' @export
topleft_legend <- function(legendmargin = margin(6, 6, 6, 6), ...){
	legend <- theme(legend.position = c(.01, .99), legend.justification = c('left', 'top'), legend.box.just = "left", legend.margin = legendmargin, ...)
	return(legend)
}

#' @name ktplots
#' @param legendmargin margin(float/int, float/int, float/int, float/int)
#' @param ... passed to ggplot2::theme
#' @export
bottomleft_legend <- function(legendmargin = margin(6, 6, 6, 6), ...){
	legend <- theme(legend.position = c(.01, .01), legend.justification = c('left', 'bottom'), legend.box.just = "left", legend.margin = legendmargin, ...)
	return(legend)
}

#' @name ktplots
#' @param legendmargin margin(float/int, float/int, float/int, float/int)
#' @param ... passed to ggplot2::theme
#' @export
bottomright_legend <- function(legendmarging = margin(6, 6, 6, 6), ...){
	legend <- theme(legend.position = c(.99, .01), legend.justification = c('right', 'bottom'), legend.box.just = "left", legend.margin = legendmargin, ...)
	return(legend)
}

#' @name ktplots
#' @param fontsize float/int
#' @param linethickness float/int
#' @param ... passed to ggplot2::theme
#' @export
small_axis <- function(fontsize=8, linethickness=0.1, ...){
	axis <- theme(text = element_text(size=fontsize), axis.text = element_text(size=fontsize), axis.line = element_line(size = linethickness), axis.ticks = element_line(size = linethickness), ...)
	return(axis)
}

