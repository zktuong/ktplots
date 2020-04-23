#' shortcut to initiate building of ktplots
#'
#' @name ktplots
#' @return for internal use only. runs init ktplots
#' @examples
#' \donttest{
#' ktplots()
#' }
#' @export
ktplots <- function()
{
	setwd("~/Documents/GitHub/ktplots")
	library(ktplots)
	library(roxygen2)
	init_ktplots()
}