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
    requireNamespace("devtools")
    requireNamespace("styler")
    styler::style_pkg(".", style = styler::tidyverse_style, indent_by = 4)
    devtools::document()
    setwd("..")
    devtools::install("ktplots", dependencies = FALSE)
    setwd("ktplots")
}

#' @export

init <- function(package, dependencies = FALSE) {
    requireNamespace("styler")
    setwd(paste0("~/Documents/GitHub/", package))
    styler::style_pkg(".", style = styler::tidyverse_style, indent_by = 4)
    requireNamespace("devtools")
    devtools::document()
    setwd("..")
    devtools::install(package, dependencies = dependencies)
    setwd(package)
}
