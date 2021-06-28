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
	requireNamespace('devtools')
	devtools::document()
	setwd('..')
	devtools::install('ktplots')
	setwd('ktplots')
}

#' @export

init <- function(package)
{   
    setwd(paste0("~/Documents/GitHub/", package))
    requireNamespace('roxygen2')
    requireNamespace('devtools')
    devtools::document()
    setwd('..')
    devtools::install(package)
    setwd(package)
}
