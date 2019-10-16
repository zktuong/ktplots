#' documenting ktplots shortcut
#'
#' @return quickly redocumenting ktplots
#' @examples
#' init_ktplots() # leaves folder, install, and change back to folder
#' @import devtools
#' @export
init_ktplots <- function() {
devtools::document()
setwd('..')
devtools::install('ktplots')
setwd('ktplots')
}
