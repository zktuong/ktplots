#' shortcut to initiate building of ktplots
#' 
#' @return runs init ktplots
#' @examples
#' ktplots()
#' @export
#'      
            
ktplots <- function() 
{
	setwd("~/Documents/GitHub/ktplots")
	library(ktplots)
	library(roxygen2)
	init_ktplots()
}