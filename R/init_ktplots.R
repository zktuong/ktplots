#' documenting ktplots shortcut
#'
#' @name ktplots
#' @return for internal use only. quickly redocumenting ktplots
#' @examples
#' \donttest{
#' init_ktplots() # leaves folder, install, and change back to folder
#' }
#' @import devtools
#' @export
init_ktplots <- function() {
	devtools::document()
	setwd('..')
	devtools::install('ktplots')
	setwd('ktplots')
}
