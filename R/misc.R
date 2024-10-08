#' miscellaneous functions
#'
#' @return miscellaneous functions

#' @name misc
#' @param x a vector of values to scale values from 0 to 1
#' @examples
#' x <- range01(runif(100))
#' @export
range01 <- function(x) {
    (x - min(x))/(max(x) - min(x))
}


#' @name misc
#' @rdname misc
#' @usage x \%nin\% y
#' @param x a vector
#' @param y a vector
#' @examples
#' \donttest{
#' data$width <- data$width %nin% params$width
#' }
#' @export
"%nin%" <- function(x, y) {
    if (!is.null(x))
        x else y
}

#' @import dplyr
#' @import ggplot2
#' @name misc
#' @param fontsize float/int
#' @param ... passed to ggplot2::theme
#' @examples
#' \donttest{
#' g + small_legend()
#' }
#' @export
small_legend <- function(fontsize = 5, keysize = 0.1, marginsize = c(-0.1, 0, 0,
    0), ...) {
    small_legend_theme <- theme(legend.title = element_text(size = fontsize), legend.text = element_text(size = fontsize),
        legend.key.size = unit(keysize, "lines"), legend.margin = margin(marginsize[1],
            marginsize[2], marginsize[3], marginsize[4], unit = "cm"), ...)
    return(small_legend_theme)
}

#' @name misc
#' @param guidesize float/int
#' @param ... passed to ggplot2::theme
#' @examples
#' \donttest{
#' g + small_guide()
#' }
#' @export
small_guide <- function(guidesize = 1, ...) {
    small_guide <- guides(shape = guide_legend(override.aes = list(size = guidesize)),
        color = guide_legend(override.aes = list(size = guidesize)), ...)
    return(small_guide)
}

#' @name misc
#' @param fontsize float/int
#' @param linethickness float/int
#' @param ... passed to ggplot2::theme
#' @examples
#' \donttest{
#' g + small_axis()
#' }
#' @export
small_axis <- function(fontsize = 4, linethickness = 0.1, ...) {
    axis <- theme(text = element_text(size = fontsize), axis.text = element_text(size = fontsize),
        axis.text.x = element_text(size = fontsize), axis.text.y = element_text(size = fontsize),
        axis.line = element_line(linewidth = linethickness), axis.ticks = element_line(linewidth = linethickness),
        ...)
    return(axis)
}

#' @name misc
#' @param linethickness float/int
#' @param ... passed to ggplot2::theme
#' @examples
#' \donttest{
#' g + small_grid()
#' }
#' @export
small_grid <- function(linethickness = 0.1, panelthickness = 0.3, ...) {
    grid <- theme(panel.grid = element_line(linewidth = linethickness), panel.border = element_rect(linewidth = panelthickness),
        ...)
    return(grid)
}

#' @name misc
#' @param legendmargin margin(float/int, float/int, float/int, float/int)
#' @param ... passed to ggplot2::theme
#' @examples
#' \donttest{
#' g + topright_legend()
#' }
#' @export
topright_legend <- function(legendmargin = margin(6, 6, 6, 6), ...) {
    legend <- theme(legend.position = c(0.99, 0.99), legend.justification = c("right",
        "top"), legend.box.just = "right", legend.margin = legendmargin, ...)
    return(legend)
}

#' @name misc
#' @param legendmargin margin(float/int, float/int, float/int, float/int)
#' @param ... passed to ggplot2::theme
#' @examples
#' \donttest{
#' g + topleft_legend()
#' }
#' @export
topleft_legend <- function(legendmargin = margin(6, 6, 6, 6), ...) {
    legend <- theme(legend.position = c(0.01, 0.99), legend.justification = c("left",
        "top"), legend.box.just = "left", legend.margin = legendmargin, ...)
    return(legend)
}

#' @name misc
#' @param legendmargin margin(float/int, float/int, float/int, float/int)
#' @param ... passed to ggplot2::theme
#' @examples
#' \donttest{
#' g + bottomleft_legend()
#' }
#' @export
bottomleft_legend <- function(legendmargin = margin(6, 6, 6, 6), ...) {
    legend <- theme(legend.position = c(0.01, 0.01), legend.justification = c("left",
        "bottom"), legend.box.just = "left", legend.margin = legendmargin, ...)
    return(legend)
}

#' @name misc
#' @param legendmargin margin(float/int, float/int, float/int, float/int)
#' @param ... passed to ggplot2::theme
#' @examples
#' \donttest{
#' g + bottomright_legend()
#' }
#' @export
bottomright_legend <- function(legendmargin = margin(6, 6, 6, 6), ...) {
    legend <- theme(legend.position = c(0.99, 0.01), legend.justification = c("right",
        "bottom"), legend.box.just = "left", legend.margin = legendmargin, ...)
    return(legend)
}